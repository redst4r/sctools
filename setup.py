from setuptools import setup, find_packages

setup(name='sctools',
      version=0.1,
      description='Tools that come handy when working in scanpy all day long',
      url='http://github.com/redst4r/sctools/',
      author='redst4r',
      maintainer='redst4r',
      maintainer_email='redst4r@web.de',
      license='GNU GPL 3',
      keywords='scanpy, scrnaseq',
      packages=find_packages(),
      install_requires=[
          'toolz',
          'numpy',
          'tqdm',
          'scipy',
          'pandas',
          'scanpy',
          'plotnine',
          'rnaseqtools @git+https://github.com/redst4r/rnaseqtools',
          'scrublet',
          'harmonypy'
          ],
      zip_safe=False)
