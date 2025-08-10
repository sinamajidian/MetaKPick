
from distutils.core import setup
setup(
  name = 'metakpick',
  packages = ['metakpick'],
  version = '0.0.1',
  license='MIT',
  description = 'The toolkits to recover the IG region from LCL dataset.',
  author = 'Sina Majidian',
  author_email = 'sina.majidian@gmail.com',
  url = 'https://github.com/sinamajidian/MetaKPick',
  download_url = 'https://github.com/sinamajidian/MetaKPick/tarball/master',
  keywords = ['metagenomics', 'classification', 'machine learning', 'random forest'],
  install_requires=[
          'numpy',
          'pandas',
          'matplotlib',
          'scikit-learn',
      ],
  include_package_data=True,
  #package_data={'metakpick': ['files/*', 'files/nodes.dmp', 'files/seqid2taxid.map_tax_uniq']},
  zip_safe = False,
  classifiers=[
    'Programming Language :: Python :: 3',
  ],
  entry_points={"console_scripts": ["metakpick = metakpick.metakpick:main"],},
)