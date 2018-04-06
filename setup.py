try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='OptAux',
      version='0.01',
      description='Code to run OptAux',
      author='Colton Lloyd',
      author_email='cjlloyd@ucsd.edu',
      url='https://github.com/coltonlloyd/OptAux',
      install_requires=[],
      packages=['OptAux'],
      )