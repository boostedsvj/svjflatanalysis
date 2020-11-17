from setuptools import setup
setup(
    name          = 'svjflatanalysis',
    version       = '0.8',
    license       = 'BSD 3-Clause License',
    description   = 'Description text',
    url           = 'https://github.com/tklijnsma/svjflatanalysis.git',
    download_url  = 'https://github.com/tklijnsma/svjflatanalysis/archive/v0_8.tar.gz',
    author        = 'Thomas Klijnsma',
    author_email  = 'tklijnsm@gmail.com',
    packages      = ['svjflatanalysis'],
    zip_safe      = False,
    scripts       = [],
    include_package_data = True,
    install_requires = ['seutils', 'uproot'],
    )