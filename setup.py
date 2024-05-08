from setuptools import Extension, setup

module = Extension("mysymnmf", sources = ['symnmfmodule.c','symnmf.c'])
setup(name = 'mysymnmf', version='1.0', description='test',ext_modules=[module],headers=['symnmf.h'])