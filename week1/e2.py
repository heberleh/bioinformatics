
import sys # you must import "sys" to read from STDIN

def patternCount(text, code):
    c = 0
    k = len(code)
    for i in range(len(text)-k+1):        
        if text[i:i+k] == code:
            c += 1
    return c

def freqWords(text, k):
    freq = {}
    max_count = 0
    for i in range(len(text)-k+1):
        word = text[i:i+k]        
        if word not in freq:            
            count = patternCount(text, word)
            freq[word] = count
            if count > max_count:
                max_count = count
    words = set()    
    for word in freq:
        if freq[word] == max_count:
            words.add(word)
    return words

lines = sys.stdin.read().splitlines() # read in the input from STDIN

code = lines[0]
k = int(lines[1])
solution = None
if len(lines) > 2:
    solution = lines[2]

result_unformated = freqWords(code,k)
result = ""
for word in result_unformated:
    result += word + " "
result = result[:-1]
print(result)