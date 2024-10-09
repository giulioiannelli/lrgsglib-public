import time
import random
import string

class Chronometer:
    chronometer_data = {}  # Class variable to store accumulated timing data

    def __init__(self, name, auto_print=False):
        self.name = name
        self.id = ''.join(random.choices(string.ascii_letters + string.digits, k=8))
        self.start_time = time.perf_counter()
        self.end_time = None
        self.auto_print = auto_print

    def end(self):
        if self.end_time is None:
            self.end_time = time.perf_counter()
            elapsed_time = self.get_elapsed_time()
            if self.name not in Chronometer.chronometer_data:
                Chronometer.chronometer_data[self.name] = {"total_time": 0, "count": 0}
            Chronometer.chronometer_data[self.name]["total_time"] += elapsed_time
            Chronometer.chronometer_data[self.name]["count"] += 1
            if self.auto_print:
                print(f"Chronometer '{self.name}' (ID: {self.id}) stopped. Elapsed time: {elapsed_time:.4g} seconds")
        else:
            if self.auto_print:
                print(f"Chronometer '{self.name}' with ID '{self.id}' has already been stopped.")

    def get_elapsed_time(self):
        if self.end_time is not None:
            elapsed_time = self.end_time - self.start_time
        else:
            elapsed_time = time.perf_counter() - self.start_time
        return elapsed_time

    @staticmethod
    def print_all_chronometers():
        print("\nSummary of all Chronometers:")
        header = f"{'Function':<25} {'Total Time (s)':<20} {'Calls':<10} {'Average Time (s)':<20}"
        print(header)
        print("-" * len(header))
        sorted_data = sorted(Chronometer.chronometer_data.items(), key=lambda item: item[1]['total_time'] / item[1]['count'], reverse=True)
        for name, data in sorted_data:
            average_time = data["total_time"] / data["count"]
            formatted_line = f"{name:<25} {data['total_time']:<20.4g} {data['count']:<10} {average_time:<20.4g}"
            print(formatted_line)