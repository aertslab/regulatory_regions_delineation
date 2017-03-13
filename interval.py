class Interval:
    @staticmethod
    def merge(intervals):
        intervals = sorted(intervals)
        if len(intervals) == 0:
            return intervals
        merged_intervals = []
        start, end = intervals[0].start, intervals[0].end
        for cur_interval in intervals[1:]:
            if cur_interval.isempty():
                continue
            cur_start, cur_end = cur_interval.start, cur_interval.end
            if cur_start <= end:
                end = max(end, cur_end)
            else:
                merged_intervals.append(Interval(start, end))
                start, end = cur_start, cur_end
        # Changed on 27/06/2011: do not add empty intervals ...
        if start < end:
            merged_intervals.append(Interval(start, end))
        return merged_intervals

    def __init__(self, start, end):
        self.start = start
        self.end = end if start < end else start

    def __str__(self):
        if self.isempty():
            return "{}"
        else:
            return "[{0:d},{1:d}[".format(self.start, self.end)

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __contains__(self, element):
        return element >= self.start and element < self.end

    def isempty(self):
        return self.start >= self.end

    def overlaps_with(self, other):
        return self.start < other.end and self.end > other.start

    def __lt__(self, other):
        if self.start < other.start:
            return True
        elif self.start > other.start:
            return False
        elif self.end < other.end:
            return True
        elif self.end > other.end:
            return False
        else:
            return False
