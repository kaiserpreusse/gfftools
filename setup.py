from setuptools import setup, find_packages

from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='gfftools',
      use_scm_version={
          "root": ".",
          "relative_to": __file__,
          "local_scheme": "node-and-timestamp"
      },
      setup_requires=['setuptools_scm'],
      description='Parse GFF files.',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='https://github.com/kaiserpreusse/gfftools',
      author='Martin Preusse',
      author_email='martin.preusse@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=[],
      keywords=['GFF', 'Genomics'],
      zip_safe=False,
      classifiers=[
          'Programming Language :: Python',
          'Intended Audience :: Developers',
      ],
      )
