from tools.entry import entrypoint
import os

if __name__ == "__main__":
    dir = os.path.dirname(os.path.realpath(__file__))
    entrypoint(dir)
