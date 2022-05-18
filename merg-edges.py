from pathlib import Path
import pandas

cnt = 0
frames = []
pathlist = Path("./merg/").glob('**/*.pkl')
for path in pathlist:
    frame = pandas.read_pickle(path)
    print(f"Records in this pickle: {len(frame)}")
    cnt += len(frame)
    print(f"Current cnt: {cnt}")
    frames.append(frame)

rez = pandas.concat(frames).sort_index()
print(f"Total records: {len(rez)}")
print(rez.head(50))

rez.to_pickle("./data/all_my_systems.pkl")