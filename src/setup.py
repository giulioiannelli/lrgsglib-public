from setuptools import setup

setup(
   name='LRGSG',
   version='0.1',
   description='A useful module',
   author='Man Foo',
   author_email='foomail@foo.example',
   packages=['LRGSG_package'],  #same as name
   install_requires=['numpy'], #external packages as dependencies
)
