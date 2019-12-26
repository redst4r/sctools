from setuptools import setup

setup(name='sctools',
      version=0.1,
      description='Tools that come handy when working in scanpy',
      url='http://github.com/redst4r/sctools/',
      author='redst4r',
      maintainer='redst4r',
      maintainer_email='redst4r@web.de',
      license='GNU GPL 3',
      keywords='scanpy, scrnaseq',
      packages=['sctools'],
      install_requires=[
          'toolz',
          'numpy',
          'tqdm',
          'scipy',
          'pandas',
          'scanpy',
          # 'rnaseqtools'  # TODO turn this into a github dependency
          ],
      zip_safe=False)
