
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import rcParams
# 和用 plt.rcParams or matplotlib.rcParams 是一样的
import sys
if sys.platform[:3] == 'win':
    dict_rc = {
        'pdf.fonttype': 42,
        'font.family': ['Arial']
    }
elif sys.platform[:5] == 'linux':
    dict_rc = {
        'pdf.fonttype': 42,
        'font.family': ['Arial']
    }
else:
    dict_rc = {
        'pdf.fonttype': 42,
    }
rcParams.update(dict_rc)

# 2024.06.03. offtable 添加 threshold 分界线，默认为 None，常用的是 2

def offtable(offtargets, target_guide,  length_pam = 3, 
                         col_seq='best_target', col_score='track_score', col_mismatch='mismatch', col_loc='target_location',
                         title=None, font='Arial', font_size=9, 
                         box_size_x=15, box_size_y=20, box_gap=1, threshold=None,
                         x_offset=15, y_offset=35, dpi=300, savefig=None):
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
        if threshold is not None:
            n_positive = sum(offtargets[col_score]>=threshold)
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
    # dpi=300
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
        # 2025.06.05. 如果有负数的，用红色显示
        if seq[col_score]>0:
            text_color = 'black'
        else:
            text_color = 'red'
        ax.text(x_offset + (len(target_guide) + 2) * box_size_x, y + box_size_y / 2, round(seq[col_score],2), ha='center', va='center', family=font, fontsize=font_size, color=text_color)
        #ax.text(x_offset + (len(target_guide) + 7) * box_size_x, y + box_size_y / 2, "Target" if seq[col_mismatch] == 0 else seq[col_mismatch], ha='center', va='center', family=font, fontsize=font_size, color='red' if seq[col_mismatch] == 0 else 'black')
        ax.text(x_offset + (len(target_guide) + 4) * box_size_x, y + box_size_y / 2, seq[col_loc], ha='left', va='center', family=font, fontsize=font_size, color=text_color)

        
    # add a vertical line to indicate the PAM
    x_line = x_offset + (len(target_guide) - length_pam) * box_size_x
    y_start = y_offset # + box_size_y / 2 
    y_end = y_start + (len(offtargets)+1) * (box_size_y + box_gap)
    ax.vlines(x=x_line, ymin=y_start, ymax=y_end, color='indianred', linestyle='--')

    # 2024.06.03. add a horizontal line to indicate the threshold
    if threshold is not None:
        thresh_x_start = x_offset
        thresh_x_end = x_offset + len(target_guide) * box_size_x
        thresh_y = y_offset + (n_positive+1) * (box_size_y + box_gap) - box_gap*0.5
        ax.hlines(y=thresh_y, xmin=thresh_x_start, xmax=thresh_x_end, color='orange', linestyle='--')


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

# summary_method: mean (default) or average, max, min, stdev, dev, coverage, cov or sum.
# number_of_bins: 700 (default) or any integer above 1
# 设置 bin_size 可以用来自动调整 number_of_bins，但是如果 **properties 里有 number_of_bins，就会被覆盖
def igv_tracking(location, file_fw, file_rv, track_name='', track_name_loc='left', 
                 fig=None, track_position = 0, track_gap = 0.2, single_height=1, bin_size=None,
                fig_scale = 0.5, aspect_ratio = 5, ax_gap = 0.02, show_title=True, spine_width=0.5,
                track_color='red', ex_length = 10000, set_ymax_fw = None, set_ymin_rv = None,
                min_ymax = None,
                savefig=None, savedpi=200, **properties):
    # only for plotting tracking-seq bw files
    import pygenometracks.tracks as pygtk
    # 一般连画时，后者都会输入 track_position
    if track_position !=0 :
        show_title = False

    if fig is None:
        fig = plt.figure(figsize=(fig_scale, fig_scale))
    
    track_height=2*single_height+track_gap*2+ax_gap
    track_position = track_position - track_height
    fw_ax = fig.add_axes([0, track_position + track_gap + single_height + ax_gap, aspect_ratio, single_height])
    rv_ax = fig.add_axes([0, track_position + track_gap , aspect_ratio, single_height])

    location = location.replace(',','')
    chrom = location.split(':')[0]
    start_region,end_region = location.split(':')[1].split('-')
    start_region = int(start_region) - ex_length
    end_region = int(end_region) + ex_length

    track_config_fw  = dict(file=file_fw)
    tk_fw = pygtk.BigWigTrack(track_config_fw)
    tk_fw.properties['color'] = track_color
    tk_fw.properties['negative_color'] = track_color
    if bin_size is not None:
        n_bins = (end_region-start_region)//bin_size
        tk_fw.properties['number_of_bins'] = n_bins
    # for properties in kwargs:
    for key, value in properties.items():
        tk_fw.properties[key] = value
    tk_fw.plot(fw_ax,chrom,start_region,end_region,)
    ymax_fw = fw_ax.get_ylim()[1]
    print('ymax_fw',ymax_fw)
    if set_ymax_fw:
        fw_ax.set_ylim(0,set_ymax_fw)
        real_ymax_fw = set_ymax_fw
    else:
        if min_ymax is not None:
            ymax_fw = max(ymax_fw,min_ymax)
        fw_ax.set_ylim(0,ymax_fw)
        real_ymax_fw = ymax_fw
    fw_ax.set_xlim(start_region,end_region)
    # hide the spine and ticks 
    for spine in fw_ax.spines.values():
        spine.set_visible(False)
    fw_ax.spines['bottom'].set_visible(True)
    fw_ax.spines['bottom'].set_linewidth(fig_scale*spine_width)
    #fw_ax.spines['bottom'].set_color(track_color)
    fw_ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    track_config_rv  = dict(file=file_rv)
    tk_rv = pygtk.BigWigTrack(track_config_rv)
    tk_rv.properties['color'] = track_color
    tk_rv.properties['negative_color'] = track_color
    if bin_size is not None:
        n_bins = (end_region-start_region)//bin_size
        tk_rv.properties['number_of_bins'] = n_bins
    # for properties in kwargs:
    for key, value in properties.items():
        tk_rv.properties[key] = value
    tk_rv.plot(rv_ax,chrom,start_region,end_region,)
    ymin_rv = rv_ax.get_ylim()[0]
    print('ymin_rv',ymin_rv)
    if set_ymin_rv:
        rv_ax.set_ylim(set_ymin_rv,0)
        real_ymin_rv = set_ymin_rv
    else:
        if min_ymax is not None:
            ymin_rv = min(ymin_rv,-min_ymax)
        rv_ax.set_ylim(ymin_rv,0)
        real_ymin_rv = ymin_rv
    rv_ax.set_xlim(start_region,end_region) # 实际上没必要，因为 sharex='col'
    # hide the spine and ticks 
    for spine in rv_ax.spines.values():
        spine.set_visible(False)
    rv_ax.spines['top'].set_visible(True)
    rv_ax.spines['top'].set_linewidth(fig_scale*spine_width)
    #rv_ax.spines['top'].set_color(track_color)
    rv_ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    # add y range on the left top
    # 如果 properties 里有 summary_method = 'sum'， 则除以 bin size
    showed_ymax = real_ymax_fw
    showed_ymin = real_ymin_rv
    if 'summary_method' in properties:
        if properties['summary_method'] == 'sum':
            showed_ymax = real_ymax_fw/bin_size
            showed_ymin = real_ymin_rv/bin_size
    fw_ax.text(start_region, real_ymax_fw, f'{showed_ymin:.1f}-{showed_ymax:.1f}', ha='left', va='top', fontsize=fig_scale*15)

    # add track name to the left or right
    x_range = end_region - start_region
    x_gap = x_range*0.02
    if track_name_loc == 'left':
        fw_ax.text(start_region-x_gap, 0, track_name, ha='right', va='center', fontsize=fig_scale*20)
    else:
        fw_ax.text(end_region+x_gap, 0, track_name, ha='left', va='center', fontsize=fig_scale*20)

    print(f'{chrom}:{start_region}-{end_region}')
    if show_title:
        region_length = round((end_region - start_region)/1000,1)
        fw_ax.set_title(f'{chrom}:{start_region}-{end_region}\n({region_length:g} kb)',loc='center',fontsize=fig_scale*20)

    if savefig is not None:
        plt.savefig(savefig, bbox_inches='tight', dpi=savedpi)

    return fig, track_position


def igv_single(location, file, fig=None, track_name='', track_name_loc='left',
                track_position = 0, track_gap = 0.2, bin_size=None,
                fig_scale = 0.5, aspect_ratio = 5, show_title=True, spine_width=0.5,
                track_color='red', ex_length = 10000, set_ymax_single = None, min_ymax=None,
                savefig=None, savedpi=200, **properties):
    # for plotting a general bw file
    import pygenometracks.tracks as pygtk

    # 一般连画时，后者都会输入 track_position
    if track_position !=0 :
        show_title = False

    if fig is None:
        fig = plt.figure(figsize=(fig_scale, fig_scale))
    
    track_height=1+track_gap*2
    track_position = track_position - track_height
    single_ax = fig.add_axes([0, track_position + track_gap, aspect_ratio, 1])

    location = location.replace(',','')
    chrom = location.split(':')[0]
    start_region,end_region = location.split(':')[1].split('-')
    start_region = int(start_region) - ex_length
    end_region = int(end_region) + ex_length

    track_config_single  = dict(file=file)
    tk_single = pygtk.BigWigTrack(track_config_single)
    tk_single.properties['color'] = track_color
    tk_single.properties['negative_color'] = track_color
    if bin_size is not None:
        n_bins = (end_region-start_region)//bin_size
        tk_single.properties['number_of_bins'] = n_bins
    # for properties in kwargs:
    for key, value in properties.items():
        tk_single.properties[key] = value
    tk_single.plot(single_ax,chrom,start_region,end_region,)
    ymax_single = single_ax.get_ylim()[1]
    print('ymax_single',ymax_single)
    if set_ymax_single:
        single_ax.set_ylim(0,set_ymax_single)
        ylim_middle = set_ymax_single/2
        real_ymax = set_ymax_single
    else:
        if min_ymax is not None:
            ymax_single = max(ymax_single,min_ymax)
        single_ax.set_ylim(0,ymax_single)
        ylim_middle = ymax_single/2
        real_ymax = ymax_single
    single_ax.set_xlim(start_region,end_region)
    # hide the spine and ticks 
    for spine in single_ax.spines.values():
        spine.set_visible(False)
    single_ax.spines['bottom'].set_visible(True)
    single_ax.spines['bottom'].set_linewidth(fig_scale*spine_width)
    #single_ax.spines['bottom'].set_color(track_color)
    single_ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

    # add y range on the left top
    # 如果 properties 里有 summary_method = 'sum'， 则除以 bin size
    showed_ymax = real_ymax
    if 'summary_method' in properties:
        if properties['summary_method'] == 'sum':
            showed_ymax = real_ymax/bin_size
    single_ax.text(start_region, real_ymax, f'0-{showed_ymax:.0f}', ha='left', va='top', fontsize=fig_scale*15)

    # add track name to the left or right
    x_range = end_region - start_region
    x_gap = x_range*0.02
    if track_name_loc == 'left':
        single_ax.text(start_region-x_gap, ylim_middle, track_name, ha='right', va='center', fontsize=fig_scale*20)
    else:
        single_ax.text(end_region+x_gap, ylim_middle, track_name, ha='left', va='center', fontsize=fig_scale*20)

    print(f'{chrom}:{start_region}-{end_region}')
    if show_title:
        region_length = round((end_region - start_region)/1000,1)
        single_ax.set_title(f'{chrom}:{start_region}-{end_region}\n({region_length:g} kb)',loc='center',fontsize=fig_scale*20)

    if savefig is not None:
        plt.savefig(savefig, bbox_inches='tight', dpi=savedpi)

    return fig, track_position


def tracking_plot(signal_L, signal_R, bin_size=100, bins=None, 
                  figsize=(10, 3), title='',
                  show_plot=True, fig=None, ax1=None, ax2=None,
                  savefig=None, save_dpi=300):
    if bins is None:
        bins=len(signal_L)
    y_max_L = signal_L[-bins:].max()
    y_mim_R = signal_R[:bins].min()

    if fig is None:
        fig = plt.figure(figsize=figsize)
        ax1 = fig.add_axes([0.0, 0.1, 0.5, 0.8])
        ax2 = fig.add_axes([0.5, 0.1, 0.5, 0.8])

    # plot left
    ax1.plot(range(bins), signal_L[-bins:], label='Original')
    #ax1.plot(range(bins), lowess_smoothed_L[-bins:, 1], label='LOWESS', color='red')
    ax1.plot([0,bins],[0,0],label='zero',color='black')
    #ax1.plot([0,bins],[signal_threshold,signal_threshold],label='threshold_left',color='orange')
    #ax1.plot([0,bins],[-signal_threshold,-signal_threshold],label='threshold_right',color='orange')
    #ax1.plot([index_L+1,index_L+1],[y_mim_R,y_max_L],label='length cutoff',color='orange')
    ax1.set_ylim(y_mim_R,y_max_L)
    ax1.set_xlim(-1,bins+1)
    ax1.set_xlabel('distance to cleavage site (kb)')
    #ax1.set_title(left_region)

    # add xticks
    xtick_gap = 10000/bin_size # 10kb
    n_xticks = int(np.ceil(bins/xtick_gap))
    xticks = np.arange(0,n_xticks+1)*xtick_gap
    xticks_label = np.arange(0,n_xticks+1)*10
    xticks_label = np.flip(xticks_label)
    ax1.set_xticks(xticks)
    _ = ax1.set_xticklabels([f'{x:g}' for x in xticks_label])
    ax1.set_ylabel('signal difference\n(coverage per 10M reads)')

    # plot right
    ax2.plot(range(bins), signal_R[:bins], label='Original')
    #ax2.plot(range(bins), lowess_smoothed_R[:bins, 1], label='LOWESS', color='red')
    ax2.plot([0,bins],[0,0],label='zero',color='black')
    #ax2.plot([0,bins],[signal_threshold,signal_threshold],label='threshold_left',color='orange')
    #ax2.plot([0,bins],[-signal_threshold,-signal_threshold],label='threshold_right',color='orange')
    #ax2.plot([index_R,index_R],[y_mim_R,y_max_L],label='length cutoff',color='orange')
    ax2.set_ylim(y_mim_R,y_max_L)
    ax2.set_xlim(-1,bins+1)
    ax2.set_xlabel('distance to cleavage site (kb)')
    #ax2.set_title(right_region)

    # add xticks
    xtick_gap = 10000/bin_size # 10kb
    n_xticks = int(np.ceil(bins/xtick_gap))
    xticks = np.arange(0,n_xticks+1)*xtick_gap
    xticks_label = np.arange(0,n_xticks+1)*10
    ax2.set_xticks(xticks)
    _ = ax2.set_xticklabels([f'{x:g}' for x in xticks_label])

    # 左右两个图紧贴
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    ax2.set_ylabel('')

    # 人造 title
    ax2.text(0, y_max_L*1.1, title, ha='center', va='center')

    if savefig is not None:
        plt.savefig(savefig, dpi=save_dpi, bbox_inches='tight')
    
    #fig.tight_layout()
    if show_plot:
        plt.show()

    return fig, ax1, ax2


