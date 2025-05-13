# import os
# import shutil
# import pytest


# @pytest.fixture
# def prep_temp_dir(tmp_path_factory):
#     temp_dir = tmp_path_factory.mktemp("test_files")
#     example_dir = os.path.join(os.path.dirname(__file__), "..", "examples")
#     for item in os.listdir(example_dir):
#         s = os.path.join(example_dir, item)
#         d = os.path.join(temp_dir, item)
#         if os.path.isfile(s):
#             shutil.copy2(s, d)
#         else:
#             shutil.copytree(s, d)
#     yield temp_dir
#     shutil.rmtree(temp_dir)
