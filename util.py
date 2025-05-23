def group_by_duplicates(iterator, equals_fn):
    result_list = []
    result_slices = []
    start_idx = 0
    first_item = next(iterator)  # Get the first item from the iterator
    result_list.append(first_item)  # Always append the first item
    prev_item = first_item  # Set prev_item to the first item
    
    for i, item in enumerate(iterator, start=1):  # Start enumeration from 1
        if not equals_fn(item, prev_item):
            result_slices.append(slice(start_idx, i))
            start_idx = i
            result_list.append(item)
            prev_item = item
    
    # Append the last slice for the final group
    result_slices.append(slice(start_idx, -1))
    
    return result_list, result_slices
    
def list_iter(ref_list, slice_obj):
    return (ref_list[index] for index in range(slice_obj.start, slice_obj.stop))
    
def find_first(iterable, predicate):
    for item in iterable:
        if predicate(item):
            return item
    return None

def pretty_format_degree(deg) -> str:
    positive = deg >= 0
    sign = '' if positive else '-'
    deg = abs(deg)
    degrees = int(deg)
    minutes = int((deg - degrees) * 60)
    seconds = ((deg - degrees) * 60 - minutes) * 60

    return f"{sign}{degrees}°{minutes}'{seconds:.5f}\""
