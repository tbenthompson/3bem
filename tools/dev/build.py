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

from tools.dev.util import oname
import os
from tools.dev.fabricate import run, after

def determine_targets(c):
    targets = dict()
    if len(c['command_params']) > 0:
        for entry in c['command_params']:
            if not entry.startswith('-'):
                targets[entry] = c['targets'][entry]
    if len(targets) > 0:
        return targets
    else:
        return c['targets']

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
    targets = determine_targets(c)
    buckets = priority_groupings(targets)
    setup_tree(c)
    for b in sorted(buckets.keys()):
        for t in buckets[b]:
            print('\nCompiling target: ' + t['binary_name'])
            compile(c, t)
    after()
    for b in sorted(buckets.keys()):
        for t in buckets[b]:
            print('\nLinking target: ' + t['binary_name'])
            link(c, t)
        after()

def setup_tree(c):
    def setup_dir(dirname):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
    setup_dir(c['build_dir'])
    for d in c['subdirs']:
        setup_dir(os.path.join(c['build_dir'], d))

def compile(c, target):
    for source in target['sources']:
        compile_runner(c['build_dir'], c['compiler'], source, target['cpp_flags'])

def compile_runner(build_dir, compiler, source, flags):
    run(
        compiler,
        '-c',
        source + '.cpp',
        '-o',
        oname(build_dir, source + '.o'),
        flags
    )

def link(c, t):
    objs = [oname(c['build_dir'], s + '.o')
        for s in t['sources'] + t['linked_sources']]
    binary_path = oname(c['build_dir'], t['binary_name'])
    link_runner(c['compiler'], binary_path, objs, t['link_flags'])

def link_runner(compiler, binary_path, objs, flags):
    run(compiler, '-o', binary_path, objs, flags)
