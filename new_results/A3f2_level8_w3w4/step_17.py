import json
base='<project>/research/08_theory/A3f_level8_relations'
t=json.load(open(base+'/transcript.json'))
for i in (64,65,66,67,68):
    e=t[i]
    c=e['content']
    if isinstance(c,list):
        for item in c:
            if 'toolResult' in item:
                for cc in item['toolResult']['content']:
                    print(f'--- entry {i} toolResult ---')
                    print(cc['text'][:3000])
            elif 'text' in item:
                print(f'--- entry {i} text ---'); print(item['text'][:2000])
            elif 'reasoningContent' in item:
                pass
