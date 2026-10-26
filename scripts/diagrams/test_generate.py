import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from drawio_kit import Diagram

def test():
    print("Testing diagram imports and canvas creation...")
    d = Diagram(100, 100)
    print("Success!")

if __name__ == "__main__":
    test()
