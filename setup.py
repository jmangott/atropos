import os
import pathlib
import subprocess

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as build_ext_orig
from setuptools.command.install import install

here = os.path.abspath(os.path.dirname(__file__))

class MyInstall(install):
    def run(self):
        try:
            cmake_args = ["-DCMAKE_BUILD_TYPE=Release"]
            subprocess.call(
                ["cmake", "-B", "build"] + cmake_args,
                cwd=os.path.join(here, "atropy/core"),
            )
            subprocess.call(
                ["cmake", "--build", "build"], cwd=os.path.join(here, "atropy/core")
            )
            print("DID WORK")
        except Exception as e:
            print(e)
            print("Error compiling")
            exit(1)
        else:
            print("DID NOT WORK")
            install.run(self)


class CMakeExtension(Extension):
    def __init__(self, name):
        # don't invoke the original build_ext for this special extension
        super().__init__(name, sources=[])


class build_ext(build_ext_orig):
    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)
        super().run()

    def build_cmake(self, ext):
        cwd = pathlib.Path().absolute()

        # these dirs will be created in build_py, so if you don't have
        # any python sources to bundle, the dirs will be missing
        build_temp = pathlib.Path(self.build_temp)
        build_temp.mkdir(parents=True, exist_ok=True)
        extdir = pathlib.Path(self.get_ext_fullpath(ext.name))
        extdir.mkdir(parents=True, exist_ok=True)

        cmake_args = [
            '-DCMAKE_BUILD_TYPE=Release'
        ]

        os.chdir(str(build_temp))
        self.spawn(['cmake', str(cwd)] + cmake_args)
        if not self.dry_run:
            self.spawn(['cmake', '--build', '.'])
        # Troubleshooting: if fail on line above then delete all possible 
        # temporary CMake files including "CMakeCache.txt" in top level dir.
        os.chdir(str(cwd))


setup(
    cmdclass={
        'install': MyInstall
    }
)
