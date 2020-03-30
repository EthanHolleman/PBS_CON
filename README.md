# PBS_CON
Code for finding primer binding site consensus sequences and visualizing. This was written in an afternoon for personal use so expect some quirks if you really try and dig into it. 

## Example Run
Get yourself a clustal format multiple sequence alignment of your sequences. This is the only required data file. A normal command might look like

```
python main.py -s 800 -e 1000 -m TGG -c ./GMR30_align
```
Where `GMR30_align` is your multible sequence alignment.

### Help Menu
```
usage: main.py [-h] [-s S] [-e E] [-m M] [-c C] [-o O]

optional arguments:
  -h, --help  show this help message and exit
  -s S        Index in sequence to begin search
  -e E        Index in sequence to end search
  -m M        Starting motif
  -c C        Path to clustal formated alignment file
  -o O        Output path for png graph. Defauly = .
  ```

## Output
- Fasta file with top PBS consensus canidates
- Png visualization of your top alignment (like the example below)

![top_canidate](https://user-images.githubusercontent.com/45807040/77862393-bc100180-71e0-11ea-82a2-e48d128381d6.png)

The colored nodes represent the most frequently occurring base at that index. The larger the node the higher the frequency, edges show possible alternative consensus sequences.

The consensus sequence that would be written to the fasta file would be **TGGTGTATGGTC**
