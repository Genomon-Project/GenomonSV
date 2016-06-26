#!/usr/bin/env python

from distutils.core import setup

setup(name='genomonsv',
      version='0.4.1',
      description='Python tools for detecting somatic structural variation from cancer genome sequencing data.',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/genomonsv',
      package_dir = {'': 'lib'},
      packages=['genomonsv'],
      scripts=['GenomonSV'],
      license='GPL-3'
     )

