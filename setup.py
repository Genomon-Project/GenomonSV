#!/usr/bin/env python

from distutils.core import setup

setup(name='GenomonSV',
      version='0.1.0',
      description='Python tools for detecting somatic structural variation from cancer genome sequencing data.',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/genomonSV',
      package_dir = {'': 'lib'},
      packages=['genomonSV'],
      scripts=['GenomonSV'],
      license='GPL-3'
     )

