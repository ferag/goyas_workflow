import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

import yaml

from scripts.final_update import update_metadata as metadata


BASE_XML = Path(__file__).resolve().parents[1] / 'base_iso19139.xml'
NS = metadata.namespaces


class TransferSizeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.data = self.directory / 'published.tif'
        self.config = {'dataset': {'file': str(self.data)}}

    def make_file(self, size):
        with self.data.open('wb') as stream:
            stream.truncate(size)

    def generate(self, stat_error=None):
        config_path = self.directory / 'config.yaml'
        config_path.write_text(yaml.safe_dump(self.config))
        session_path = self.directory / 'session.json'
        session_path.write_text(json.dumps({'cookies': {}}))
        output = self.directory / 'metadata.xml'
        response = {
            'tif_response': {'url': 'https://storage.example/published.tif',
                             'filename': 'published.tif'},
            'png_response': {'url': 'https://catalog.example/thumbnail.png'},
        }
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.ExitStack() as stack:
            stack.enter_context(contextlib.redirect_stdout(stdout))
            stack.enter_context(contextlib.redirect_stderr(stderr))
            session = stack.enter_context(patch.object(metadata.requests, 'Session'))
            session.return_value.put.return_value.status_code = 200
            if stat_error is not None:
                stack.enter_context(patch.object(metadata.os.path, 'getsize', side_effect=stat_error))
            metadata.update_metadata(
                BASE_XML, 'test-record', response, output,
                'https://maps.example/wms', config_path, session_path,
            )
        root = ET.parse(output).getroot()
        # The remote update must receive the same metadata as the local output.
        remote = ET.fromstring(session.return_value.put.call_args.kwargs['data'])
        self.assertEqual(ET.tostring(root), ET.tostring(remote))
        return root, stdout.getvalue(), stderr.getvalue()

    def test_mib_conversion_and_two_decimal_places(self):
        for size, expected in [(1048576, '1.00'), (18876404, '18.00'),
                               (1294991, '1.23'), (1300000, '1.24'), (0, '0.00')]:
            with self.subTest(size=size):
                self.make_file(size)
                with contextlib.redirect_stdout(io.StringIO()) as log:
                    self.assertEqual(metadata.get_data_file_size_mb(self.config), expected)
                self.assertIn(f'Data file size: {expected} MiB ({size} bytes)', log.getvalue())

    def test_size_is_in_download_transfer_before_online(self):
        self.make_file(18876404)
        root, log, warnings = self.generate()
        self.assertEqual(warnings, '')
        self.assertIn('Data file size: 18.00 MiB (18876404 bytes)', log)
        sizes = root.findall('.//gmd:transferSize', NS)
        self.assertEqual(len(sizes), 1)
        transfers = root.findall('.//gmd:MD_DigitalTransferOptions', NS)
        self.assertEqual(len(transfers), 2)
        for transfer in transfers:
            url = transfer.findtext('gmd:onLine/gmd:CI_OnlineResource/gmd:linkage/gmd:URL', namespaces=NS)
            size = transfer.find('gmd:transferSize', NS)
            if url == 'https://storage.example/published.tif':
                self.assertIs(size, sizes[0])
                self.assertEqual(size.findtext('gco:Real', namespaces=NS), '18.00')
                self.assertEqual(list(transfer)[0], size)
                self.assertEqual(list(transfer)[1].tag, f"{{{NS['gmd']}}}onLine")
            else:
                self.assertEqual(url, 'https://maps.example/wms')
                self.assertIsNone(size)
        self.assertIsNone(root.find('.//gmd:MD_BrowseGraphic/gmd:transferSize', NS))

    def test_missing_file_warns_and_still_generates_and_uploads(self):
        root, _, warnings = self.generate()
        self.assertIn('Warning: Cannot determine data file size', warnings)
        self.assertEqual(root.findall('.//gmd:transferSize', NS), [])
        self.assertIn('https://storage.example/published.tif',
                      [node.text for node in root.findall('.//gmd:URL', NS)])

    def test_stat_permission_failure_is_nonfatal(self):
        root, _, warnings = self.generate(PermissionError('permission denied'))
        self.assertIn('permission denied', warnings)
        self.assertEqual(root.findall('.//gmd:transferSize', NS), [])

    def test_missing_path_warns_and_legacy_path_is_supported(self):
        with contextlib.redirect_stderr(io.StringIO()) as warnings:
            self.assertIsNone(metadata.get_data_file_size_mb({}))
        self.assertIn('Warning:', warnings.getvalue())
        self.make_file(1048576)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(metadata.get_data_file_size_mb({'file': str(self.data)}), '1.00')


if __name__ == '__main__':
    unittest.main()
