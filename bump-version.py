#!/usr/bin/env python
# 
#  bump-version.py
#  SEDMachine
#  
#  Created by Alexander Rudy on 2011-10-04.
#  Copyright 2011 Alexander Rudy. All rights reserved.
# 

import sys,os,re,shutil
from tempfile import mkstemp

if __name__ == '__main__':
    
    VERSION_FILENAME = "VERSION"
    
    VERSION_FILE = open(VERSION_FILENAME)
    
    VERSION = VERSION_FILE.read()
    
    VERSION_PARTS = VERSION.split(".")
    
    VERSION_RE = re.compile("(#\s+Version:\s+)[\d\.]+")
    
    if len(sys.argv) < 2:
        
        lastNumber = int(VERSION_PARTS[-1])
                
        VERSION_PARTS[-1] = lastNumber + 1
                
        VERSION = ""
        
        for part in VERSION_PARTS:
            VERSION += str(part) + "."
            
        VERSION = VERSION.rstrip(".")
        
    else:
        
        VERSION = sys.argv[1]
        
    print("ESTABLISHING VESION %s" % VERSION)
    
    for root,dirs,files in os.walk("."):
        if ".git" in dirs:
            dirs.remove(".git")
        for afile in files:
            if afile != afile.rstrip(".py") and afile != os.path.basename(__file__):
                fh, abs_path = mkstemp()
                new_file = open(abs_path,"w")
                old_file = open(root+"/"+afile)
                for line in old_file:
                    result = VERSION_RE.search(line)
                    if result:
                        output = result.groups()[0] + VERSION + "\n"
                    else:
                        output = line
                    new_file.write(output)
                new_file.close()
                os.close(fh)
                old_file.close()
                os.remove(root+"/"+afile)
                shutil.move(abs_path,root+"/"+afile)
    