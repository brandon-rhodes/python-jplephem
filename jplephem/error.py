#!/usr/bin/env python

class InvalidTimestampError(ValueError):
    def __init__(self, message, min_timestamp, max_timestamp):
        self.message = message
        self.min_timestamp = min_timestamp
        self.max_timestamp = max_timestamp
