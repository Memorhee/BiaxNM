@echo off
cd /d "%~dp0"

python src\main.py examples\rect_aci318-25.toml
python src\main.py examples\circle_kds14-22.toml
python src\main.py examples\hollow_kds14-22-app.toml
python src\main.py examples\track_ec2-05-2.toml
python src\main.py examples\phc_kds24-25.toml
python src\main.py examples\src_ec2-04-1-1.toml

pause
