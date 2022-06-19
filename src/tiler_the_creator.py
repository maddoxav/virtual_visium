OPENSLIDE_PATH = r'C:\Users\maddoxav\openslide-win64-20171122\openslide-win64-20171122\bin'
import os
with os.add_dll_directory(OPENSLIDE_PATH):
    import openslide 