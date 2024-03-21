import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

def offtable(offtargets, target_guide, 
                         col_seq='best_target', col_score='track_score', col_mismatch='mismatch', col_loc='target_location',
                         title=None, font='Arial', font_size=9, 
                         box_size_x=15, box_size_y=20, box_gap=1,
                         x_offset=15, y_offset=35, dpi=100, savefig=None):
    # Facecolor
    color_dict = {
        'A': 'lightgreen',
        'T': 'lightblue',
        'C': 'lightcoral',
        'G': 'lightgoldenrodyellow',
        'N': 'lightgrey',
        '—': 'orange',
        '-': 'orange'
    }

    # If offtargets is a DataFrame, convert to list of dictionaries
    if isinstance(offtargets, pd.DataFrame):
        offtargets = offtargets.to_dict(orient='records')

    # Configuration
    # title=None
    # font='Arial'
    # font_size = 9
    # box_size_x = 15 # 一个碱基图形的宽度
    # box_size_y = 20 # 一个碱基图形的高度
    # box_gap = 1 # 两行之间的间隔
    # x_offset = 15
    # y_offset = 35
    # dpi=100
    # col_seq='best_target'
    # col_score='track_score'
    # col_mismatch='mismatch'
    # col_loc='target_location'
    width = box_size_x * (len(target_guide) + 15)
    height = y_offset + (len(offtargets) + 2) * (box_size_y + box_gap)
    fig = plt.figure(figsize=(width / 100.0, height / 100.0), dpi=dpi)
    ax = fig.add_subplot(111)

    # Plot a title
    ax.text(x_offset, 25, "Off-targets table" if title is None else f"{title}", fontsize=14, family=font)

    # Plot the reference sequence
    for i, c in enumerate(target_guide):
        x = x_offset + i * box_size_x
        y = y_offset
        base_color = color_dict.get(c, 'purple') # Default to purple if base is not recognized
        ax.add_patch(patches.Rectangle((x, y), box_size_x, box_size_y, facecolor=base_color))
        ax.text(x + box_size_x / 2, y + box_size_y / 2, c, ha='center', va='center', family=font, fontsize=font_size)
    # add column annotations
    ax.text(x_offset + (len(target_guide) + 2) * box_size_x, y_offset + box_size_y / 4, 'Track\nScore', ha='center', va='center', family=font, fontsize=font_size*1.1)
    #ax.text(x_offset + (len(target_guide) + 7) * box_size_x, y_offset + box_size_y / 2, 'Mismatch', ha='center', va='center', family=font, fontsize=font_size*1.1)
    ax.text(x_offset + (len(target_guide) + 4) * box_size_x, y_offset + box_size_y / 2, 'Coordinates', ha='left', va='center', family=font, fontsize=font_size*1.1)

    # Plot aligned sequences
    # 目前有个bug：脱靶序列如果有 insertion，长度会不一致，而且也没想到画图怎么画，只能是默认删掉第一个碱基
    for j, seq in enumerate(offtargets):
        y = y_offset + (j + 1) * (box_size_y + box_gap)
        # 长度不一致的情况
        len_out = len(seq[col_seq]) - len(target_guide)
        if len_out > 0:
            if len_out > 1:
                print(f"Warning: {seq[col_seq]} is {len_out} longer than {target_guide}")
            # 通过比较删除开头的碱基和最后的碱基，看哪个更接近target_guide
            delete_first = seq[col_seq][len_out:]
            delete_last = seq[col_seq][:-len_out]
            # 计算两个序列和target_guide的hamming distance
            hamming_first = sum([1 for i, c in enumerate(delete_first) if c != target_guide[i]])
            hamming_last = sum([1 for i, c in enumerate(delete_last) if c != target_guide[i]])
            # 选择hamming distance小的那个序列
            if hamming_first < hamming_last:
                seq[col_seq] = delete_first
            else:
                seq[col_seq] = delete_last
        elif len_out < 0:
            print(f"Warning: {seq[col_seq]} is {-len_out} shorter than {target_guide}")
            
        for i, c in enumerate(seq[col_seq]):
            # gap 的 - (minus sign) 太短了，所以替换成 — (em dash)
            if c == '-':
                c = '—'
            x = x_offset + i * box_size_x
            base_color = color_dict.get(c, 'purple') # Default to purple if base is not recognized
            if c == target_guide[i]:
                ax.add_patch(patches.Rectangle((x, y), box_size_x, box_size_y, facecolor='white')) # same
            elif target_guide[i] == 'N':
                ax.add_patch(patches.Rectangle((x, y), box_size_x, box_size_y, facecolor='white')) # N in target
            else:
                ax.add_patch(patches.Rectangle((x, y), box_size_x, box_size_y, facecolor=base_color))
            ax.text(x + box_size_x / 2, y + box_size_y / 2, "." if c == target_guide[i] else c, ha='center', va='center', family=font, fontsize=font_size, weight='bold')

        # Annotations for score, mismatches, and location coordinates
        ax.text(x_offset + (len(target_guide) + 2) * box_size_x, y + box_size_y / 2, round(seq[col_score],2), ha='center', va='center', family=font, fontsize=font_size)
        #ax.text(x_offset + (len(target_guide) + 7) * box_size_x, y + box_size_y / 2, "Target" if seq[col_mismatch] == 0 else seq[col_mismatch], ha='center', va='center', family=font, fontsize=font_size, color='red' if seq[col_mismatch] == 0 else 'black')
        ax.text(x_offset + (len(target_guide) + 4) * box_size_x, y + box_size_y / 2, seq[col_loc], ha='left', va='center', family=font, fontsize=font_size)

    # add a vertical line to indicate the PAM
    x_line = x_offset + (len(target_guide) - 3) * box_size_x
    y_start = y_offset # + box_size_y / 2 
    y_end = y_start + (len(offtargets)+1) * (box_size_y + box_gap)
    ax.vlines(x=x_line, ymin=y_start, ymax=y_end, color='indianred', linestyle='--')

    # Styling and save
    ax.set_xlim(0, width*1.1) # location 的文字太长了，所以要加长一点
    ax.set_ylim(height, 0)
    ax.axis('off')

    # # This will make the subplot(s) expand to fill the entire figure area, with no padding on any side. 
    # # In brief, make the plot bigger (not influence the font size)
    plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
    if savefig is not None:
        plt.savefig(savefig, dpi=dpi)
    plt.show()
    return ax





