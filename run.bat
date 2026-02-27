REM Comprobar argumento
if "%1"=="" (
	echo Use: run.bat ruta\al\config.yaml
	exit /b 1
)

REM Copiar el config proporcionado a config.yaml en la raíz
copy /Y "%1" config.yaml

REM Ejecutar Snakemake con los outputs finales como targets
snakemake "logs/run_all.done" --cores 1