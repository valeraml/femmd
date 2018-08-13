from distutils.core import setup, Extension
import numpy
import os

print(os.getcwd())
src_files = os.listdir("src/cppsrc")
source_list1 = []
for file in src_files:
	if file.endswith(".cpp"):
		source_list1.append("src/cppsrc/" + file)
print(source_list1)

source_list = ['../src/bonds.cpp', '../src/driver.cpp', '../src/Elements.cpp', '../src/integrator.cpp', '../src/interactions.cpp', '../src/interface.cpp', '../src/md3dsystem.cpp', '../src/nn.cpp',  '../src/particles.cpp', '../src/properties.cpp']

module1 = Extension('femmd', 
                sources = source_list1,
                include_dirs=[numpy.get_include()],
                extra_compile_args=["-std=c++14"])

setup (name = 'femmd',
        version = '1.0',
        description = 'femmd package',
        ext_modules = [module1])

