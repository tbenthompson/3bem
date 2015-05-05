"""
This is the build script for 3bem. It uses the fabricate.py build tool.
Think "make", but way better. It automatically determines dependencies
between build steps so that if one file changes, only the minimal number of
other files are recompiled or relinked.

Also, it's nice to write a build script in python instead of some arcane
domain specific language like "make".

More info:
http://code.google.com/p/fabricate/

This build script separates the different executables and libraries to be
constructed into "targets", which are compiled individually.

Some basic command line parsing is done to determine which target to compile
and whether to do a debug, release, etc build. At some point, the
sophistication of the command line parsing may need to improve.
"""
from __future__ import print_function

import os
import shutil
import sys
import copy
from build_tools.fabricate import run, after

def oname(build_dir, filename):
    return os.path.join(build_dir, filename)

def build_target(c, target):
    setup_tree(c)
    compile(c, target)
    after()
    link(c, target)
    after()

def setup_tree(c):
    def setup_dir(dirname):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    setup_dir(c['build_dir'])
    for k in c['subdirs']:
        setup_dir(os.path.join(c['build_dir'], c['subdirs'][k]))

def compile(c, target):
    to_compile = copy.copy(target['sources'])
    for source in to_compile:
        compile_runner(c['build_dir'], c['compiler'], source,
                       target['cpp_flags'])

def compile_runner(build_dir, compiler, source, flags):
    args = [
        compiler,
        '-c',
        source + '.cpp',
        '-o',
        oname(build_dir, source + '.o'),
    ]
    args.extend(flags)
    run(*args)

def link(c, t):
    obj_raw_names = t['sources'] + t['linked_sources']
    objs = [oname(c['build_dir'], s + '.o') for s in obj_raw_names]
    binary_path = oname(c['build_dir'], t['binary_name'])
    link_runner(c['compiler'], binary_path, objs, t['link_flags'])

def link_runner(compiler, binary_path, objs, flags):
    run(compiler, '-o', binary_path, objs, flags)
