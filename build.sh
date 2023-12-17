rm -fr build dist emprove.egg-info venv
/usr/bin/python3.8 -m venv venv
source venv/bin/activate
/usr/bin/python3.8 setup.py bdist_wheel
/usr/bin/python3.8 -m pip install ./dist/emprove-0.1.0-cp38-cp38-linux_x86_64.whl
/usr/bin/python3.8 -m pip install ./dist/emprove-0.1.0-cp38-cp38-linux_x86_64.whl   --force-reinstall
rm ~/software/emprove

