from interval import Interval


class PieceWiseLocation:
    @staticmethod
    def singleton(chromosome, on_positive_strand, start, end):
        return PieceWiseLocation(chromosome, on_positive_strand, [Interval(start, end)])

    def __init__(self, chromosome, on_positive_strand, intervals):
        self.chromosome = chromosome
        self.on_positive_strand = on_positive_strand
        self.intervals = Interval.merge(intervals)

    def __str__(self):
        if self.on_positive_strand:
            string = "{0:s}: 5' {1:s}".format(self.chromosome, str(self.intervals[0]))
            for interval in self.intervals[1:]:
                string += "," + str(interval)
            string += " 3' (+, 0-based)"
            return string
        else:
            string = "{0:s}: 3' {1:s}".format(self.chromosome, str(self.intervals[0]))
            for interval in self.intervals[1:]:
                string += "," + str(interval)
            string += " 5' (-, 0-based)"
            return string

    def span(self):
        return Interval(self.intervals[0].start, self.intervals[-1].end)

    def __add__(self, other):
        return PieceWiseLocation(self.chromosome, self.on_positive_strand,
                                 self.intervals + other.intervals)

    def __sub__(self, other):
        result_intervals = self.intervals[:]
        for subtract_interval in other.intervals:
            for idx, element in enumerate(result_intervals[:]):
                if element.overlaps_with(subtract_interval):
                    result_intervals.remove(element)
                    if subtract_interval.start in element:
                        s = Interval(element.start, subtract_interval.start)
                        if not s.isempty():
                            result_intervals.append(s)
                    if subtract_interval.end in element:
                        s = Interval(subtract_interval.end, element.end)
                        if not s.isempty():
                            result_intervals.append(s)
        return PieceWiseLocation(self.chromosome, self.on_positive_strand, result_intervals)

    def isempty(self):
        # Changed on 27/06/2011: Empty intervals can also be represented as an empty list ...
        # Changed on 02/09/2011: If all intervals are empty the piecewise interval is also empty ...
        return (len(self) == 0) or all(interval.isempty() for interval in self)

    def __len__(self):
        return len(self.intervals)

    def __getitem__(self, index):
        return self.intervals[index]

    def location(self, index):
        return PieceWiseLocation(self.chromosome, self.on_positive_strand, [self.intervals[index]])

    def __iter__(self):
        for interval in self.intervals:
            yield interval

    def __contains__(self, element):
        for interval in self:
            if element in interval:
                return True
        return False

    def filter(self, element):
        return PieceWiseLocation(self.chromosome, self.on_positive_strand,
                                 list(filter(lambda interval: element in interval, self.intervals)))

    def extend_upstream(self, number_of_bases, chromosome2length):
        if number_of_bases <= 0:
            raise ValueError("Number of bases must be positive integer.")
        if self.on_positive_strand:
            first_upstream_element = self.intervals[0]
            intervals = [Interval(max(0, first_upstream_element.start - number_of_bases), first_upstream_element.end)]
            intervals.extend(self.intervals[1:])
            return PieceWiseLocation(self.chromosome, self.on_positive_strand, intervals)
        else:
            first_upstream_element = self.intervals[-1]
            intervals = self.intervals[0:-1]
            chromosome_length = chromosome2length[self.chromosome]
            intervals.append(Interval(first_upstream_element.start,
                                      min(first_upstream_element.end + number_of_bases, chromosome_length)))
            return PieceWiseLocation(self.chromosome, self.on_positive_strand, intervals)

    def extend_downstream(self, number_of_bases, chromosome2length):
        if number_of_bases <= 0:
            raise ValueError("Number of bases must be positive integer.")
        if self.on_positive_strand:
            first_downstream_element = self.intervals[-1]
            intervals = self.intervals[0:-1]
            chromosome_length = chromosome2length[self.chromosome]
            intervals.append(Interval(first_downstream_element.start,
                                      min(first_downstream_element.end + number_of_bases, chromosome_length)))
            return PieceWiseLocation(self.chromosome, self.on_positive_strand, intervals)
        else:
            first_downstream_element = self.intervals[0]
            intervals = [Interval(max(0, first_downstream_element.start - number_of_bases),
                                  first_downstream_element.end)]
            intervals.extend(self.intervals[1:])
            return PieceWiseLocation(self.chromosome, self.on_positive_strand, intervals)

    def interval_limit(self, limit):
        intervals = []
        for interval in self.intervals:
            limited_interval = Interval(max(limit.start, interval.start), min(limit.end, interval.end))
            if not limited_interval.isempty():
                intervals.append(limited_interval)
        return PieceWiseLocation(self.chromosome, self.on_positive_strand, intervals)
