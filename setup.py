try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='OptAux',
      version='0.0.1',
      description='Code to run OptAux',
      author='Colton Lloyd',
      author_email='cjlloyd@ucsd.edu',
      url='https://github.com/coltonlloyd/OptAux',
      install_requires=['openpyxl', 'pandas'],
      packages=['optaux'],
      )
