# -*- coding: utf-8 -*-
# 
#  make_directory.py
#  SEDM-Dev
#  
#  Created by Alexander Rudy on 2012-03-28.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

import os

DIRECTORIES = ["Caches","Data","Images","Logs","Partials"]

for DIR in DIRECTORIES:
    try:
        os.mkdir(DIR)
    except os.error, e:
        print "Directory %s already exists" % DIR
        