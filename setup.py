from setuptools import setup, find_packages

setup(name='sctools',
      version=0.2,
      description='Tools that come handy when working in scanpy all day long',
      url='http://github.com/redst4r/sctools/',
      author='redst4r',
      maintainer='redst4r',
      maintainer_email='redst4r@web.de',
      license='GNU GPL 3',
      keywords='scanpy, scrnaseq',
      packages=find_packages(),
      include_package_data=True,
      install_requires=[
          'toolz',
          'numpy',
          'tqdm',
          'scipy',
          'pandas',
          'scanpy>1.7',
          'plotnine',
          'rnaseqtools @git+https://github.com/redst4r/rnaseqtools',
          'scrublet',
          'harmonypy ==0.0.5',
          'gprofiler-official'
          ],
      zip_safe=False)
