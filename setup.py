from setuptools import setup

setup(
   name='subsystem-surface-code-demo',
   version='0.1.0',
   author='Asmae Benhemou',
   author_email='asmabenh@amazon.com',
   description='A basic demo to simulate the subsystem surface code',
   long_description=open('README.md').read(),
   install_requires=[
       "stim",
       "pymatching",
       "numpy",
       "sinter",
   ],
)