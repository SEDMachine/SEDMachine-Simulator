# -*- coding: utf-8 -*-
# 
#  make_directory.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2012-03-28.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

import os, shutil

from pkg_resources import resource_filename


DIRECTORIES = ["Caches","Data","Images","Logs","Partials"]

shutil.copytree(resource_filename(__name__,"Data"),"./Data",ignore=shutil.ignore_patterns("*.pyc","*.py"))

for DIR in DIRECTORIES:
    try:
        os.mkdir(DIR)
    except os.error, e:
        print "Directory %s already exists" % DIR
        