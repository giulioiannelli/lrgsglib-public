import argparse

def parse_arguments(parser: argparse.ArgumentParser) -> argparse.Namespace:
    args = parser.parse_args()
    return args

def parsDict_get(parsDict, var, key):
    return parsDict[var].get(key, None)