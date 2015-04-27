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

from tools.util import oname
import os
import shutil
import sys
import copy
from tools.fabricate import run, after

def determine_removed_targets(c):
    remove_targets = []
    for entry in c['command_params']:
        if entry.startswith('--'):
            name = entry[2:]
            remove_targets.append(name)
    return remove_targets

def determine_targets(c):
    remove_targets = determine_removed_targets(c)
    target = None
    for entry in c['command_params']:
        if entry.startswith('+'):
            if target is not None:
                print("Specify only one target or build all targets.")
                sys.exit()
            target = entry[1:]

    if target is not None:
        return dict(target = c['targets'][target]), True
    else:
        all_targets = c['targets']
        for rem in remove_targets:
            del all_targets[rem]
        return all_targets, False

def priority_groupings(targets):
    buckets = dict()
    for t in targets.values():
        p = t['priority']
        if p in buckets:
            buckets[p].append(t)
        else:
            buckets[p] = [t]
    return buckets

def run_build(c):
    targets, one_target = determine_targets(c)
    buckets = priority_groupings(targets)
    setup_tree(c)
    for b in sorted(buckets.keys()):
        for t in buckets[b]:
            c['printer']('\nCompiling target: ' + t['binary_name'])
            compile(c, t, one_target)
    after()
    for b in sorted(buckets.keys()):
        for t in buckets[b]:
            c['printer']('\nLinking target: ' + t['binary_name'])
            link(c, t)
        after()
    if 'python_wrapper' in targets:
        copy_python_wrapper_to_py(c)

def copy_python_wrapper_to_py(c):
    pylib_file = c['targets']['python_wrapper']['binary_name']
    pylib_path = oname(c['build_dir'], pylib_file)
    pylib_dest = os.path.join('py', pylib_file)
    shutil.copy(pylib_path, pylib_dest)

def setup_tree(c):
    def setup_dir(dirname):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    setup_dir(c['build_dir'])
    for k in c['subdirs']:
        setup_dir(os.path.join(c['build_dir'], c['subdirs'][k]))

def compile(c, target, one_target):
    to_compile = copy.copy(target['sources'])
    if one_target:
        to_compile.extend(target['linked_sources'])
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
