from utils import re, np

float_pattern = r'[+-]?\d*\.?\d*[eE]?[+-]?\d*'
int_pattern = r'[+-]?\d+'
delimiter = r'[\s\n]+'


def get_float(input_string, keyword):
    result = re.search(keyword + delimiter + r'(?P<value>' + float_pattern + r')', input_string)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword} in \n{input_string}')

    else:
        return float(result.group('value'))


def get_int(input_string, keyword):
    result = re.search(keyword + delimiter + r'(?P<value>' + int_pattern + r')', input_string)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword} in \n{input_string}')

    else:
        return int(result.group('value'))


def get_string(input_string, keyword):
    result = re.search(keyword + delimiter + r'(?P<value>' + r'\w+' + r')', input_string)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword} in \n{input_string}')

    else:
        return result.group('value')


def get_dict(input_string, keyword):
    result = re.search(keyword + delimiter + r'\{.*?\}', input_string, flags=re.S)
    if result is None:
        raise RuntimeError(f'keyword not found: {keyword} in \n{input_string}')
    else:
        return result.group()

def get_vector_array(input_string, is_int=False) -> np.ndarray:
    if is_int:
        pattern = re.compile(
            r'\(' + r'[\s\n]*' + int_pattern + delimiter + int_pattern + delimiter + int_pattern + r'[\s\n]*\)',
            flags=re.S)
    else:
        pattern = re.compile(
            r'\(' + r'[\s\n]*' + float_pattern + delimiter + float_pattern + delimiter + float_pattern + r'[\s\n]*\)',
            flags=re.S)
    data_list = re.findall(pattern, input_string)
    data = np.empty(shape=(len(data_list), 3))
    pattern = re.compile(r'\((.*)\)')

    if is_int:
        for i, elem in enumerate(data_list):
            data[i, :] = np.array(re.search(pattern, elem).group(1).strip().split(), dtype=int)
    else:
        for i, elem in enumerate(data_list):
            data[i, :] = np.array(re.search(pattern, elem).group(1).strip().split(), dtype=float)

    return data