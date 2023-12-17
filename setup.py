from setuptools import find_packages, setup, Extension
module1 = Extension('emprove_core', sources = ['src_cpp/emprove_core.cpp'], language='c++')

setup(
    name='emprove',
    packages=(find_packages(include=['emprove'])),
    version='0.1.0',
    description='emprove, library for assessing EM reconstruction from particles scores',
    author='Mauro Maiorca',
    license='MIT',
    ext_modules = [module1],
    package_data={
        'emprove': ['config.json'],  # Add this line to include config.json file in the emprove package
    },
    entry_points={
      'console_scripts': [
          'emprove = emprove.emprove_cmd_caller:main',
          'emprove_utils = emprove.emprove_cmd_utils:main',
          'emprove_config = emprove.emprove_config:main',
          'emprove_optimizer = emprove.emprove_cmd_optimizer:main',
        ],
      },
)

