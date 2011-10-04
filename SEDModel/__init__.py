# 
#  __init__.py
#  SEDMachine
#  
#  Created by Alexander Rudy on 2011-10-04.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 

import logging
import time
import sys

logging.basicConfig(filename="Logs/"+__name__+"-"+time.strftime("%Y-%m-%d")+".log",format="%(asctime)s:[%(levelname)-8s]:[%(name)s]:%(message)s",datefmt="%Y-%m-%d-%H:%M:%S",level=logging.DEBUG,filemode='w')