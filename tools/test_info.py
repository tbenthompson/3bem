from tools.util import files_in_dir
def get_unit_test_info(c):
    unit_test_info = dict()
    unit_test_info['main'] = dict()
    unit_test_info['main']['src'] = 'test_main'
    tests = [
        t for t in files_in_dir('test', 'cpp')
        if t.startswith('test/test_') and 'test_main' not in t
    ]
    unit_test_info['main']['lib_srcs'] = tests
    unit_test_info['main']['link_lib'] = True

    return unit_test_info

