import glob
import sys
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


def make_plot(file_name):
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['mathtext.default'] = 'regular'

  
    df = pd.read_csv(f"{file_name}")
    df = df.rename(columns=lambda x: x.strip())

    #don't allow redundant single sequence calculations to appear
    df['single'] = df['file'].str.contains('single', regex=False)
    mask = df['single']
    df = pd.concat([
        df.loc[mask].drop_duplicates(subset = ['single', df.columns[1]], keep='first'),
        df.loc[~mask], ])

    df = df.sort_values(by=['bfactor'], ascending=True)
    _ = df.plot(kind='scatter', x=df.columns[1], y=df.columns[2], c="bfactor", vmin=0.5, vmax=0.9)
    _.set_xlim([0.0,1.0])
    _.set_ylim([0.0,1.0])
    plt.savefig("{:}_scatter.png".format(file_name.split("/")[0]), dpi=300)
    plt.clf()
    plt.close()

    hist = df.hist(column=["bfactor"], bins=31, grid=False, range=[0.0,1.0])
    plt.savefig("{:}_hist.png".format(file_name.split("/")[0]), dpi=300)
    plt.clf()
    plt.close()

def combo_plot(files):
    fig, axes = plt.subplots(len(files), 1, figsize=(8,6))
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['mathtext.default'] = 'regular'

    alt_data = []
    alt_b = []
    dom_data = []
    dom_b = []
    for idx, file in enumerate(files):
        print(file)
        df = pd.read_csv(f"{file}")
        df = df.rename(columns=lambda x: x.strip())
        #don't allow redundant single sequence calculations to appear
        df['single'] = df['file'].str.contains('single', regex=False)
        mask = df['single']
        df = pd.concat([
            df.loc[mask].drop_duplicates(subset = ['single', df.columns[1]], keep='first'),
            df.loc[~mask], ])

        #create image without sufamily predictions
        # df = df[df['file'].str.contains("subfamilies")==False]
        # print("no subfamilies")
        # print(df.shape)

        axes[idx].spines['top'].set_visible(False)
        axes[idx].spines['right'].set_visible(False)

        if idx != len(files) -1:
            axes[idx].spines['bottom'].set_visible(False)
            tick_positions = [-1, -0.5, 0.0, 0.5, 1.0]
            axes[idx].set_xticklabels([])
            axes[idx].set_xticks(tick_positions)
        else:
            tick_positions = [-1, -0.5, 0.0, 0.5, 1.0]
            axes[idx].set_xticks(tick_positions)
            tick_labels = ['1.0', '0.5', '0.0', '0.5', '1.0']
            axes[idx].set_xticklabels(tick_labels)
            # pass
        axes[idx].spines['left'].set_visible(False)
        axes[idx].get_yaxis().set_ticks([])
        
        #Rfah and lymphotactin have different criteria, note this in figure
        if file.find('RFAH') == -1 and file.find('lymphotactin'):
            number_line_dict = {}
            # first column will go to the left
            name = df[df.columns[0]].to_list()
            tm_1 = df[df.columns[1]].to_list()
            tm_2 = df[df.columns[2]].to_list()
            bfactor = df["bfactor"].to_list()
            x_val = []
            b_val = []
            alt = False
            for d, a, b in zip(tm_1, tm_2, bfactor):
                if d >= a:
                    x_val.append(d)
                    b_val.append(b)
                else:
                    x_val.append(a * -1)
                    b_val.append(b)

            best_b = max(bfactor)
            
            #base the dominant fold predicted by best plddt score structure
            if tm_1[bfactor.index(best_b)] < tm_2[bfactor.index(best_b)]:
                x_val = [-1* val for val in x_val]
                alt = True

            for x, b in zip(x_val, b_val):
                if x < 0.0:
                    alt_data.append(-1*x)
                    alt_b.append(b)
                else:
                    dom_data.append(x)
                    dom_b.append(b)
            dom_name = df.columns[2].split('_')[0]
            alt_name = df.columns[1].split('_')[0]

        else:
            number_line_dict = {}
            # first column will go to the left
            name = df[df.columns[0]].to_list()
            tm_1_sub = df[df.columns[1]].to_list()
            tm_1 = df[df.columns[2]].to_list()
            tm_2_sub = df[df.columns[3]].to_list()
            tm_2 = df[df.columns[4]].to_list()

            bfactor = df["bfactor"].to_list()
            x_val = []
            b_val = []
            alt = False
            for d, a, b, d_sub, a_sub in zip(tm_1, tm_2, bfactor, tm_1_sub, tm_2_sub):
                if d_sub >= a_sub:
                    x_val.append(d)
                    b_val.append(b)
                else:
                    x_val.append(a * -1)
                    b_val.append(b)

            best_b = max(bfactor)
            
            #base the dominant fold predicted by best plddt score structure
            if tm_1[bfactor.index(best_b)] < tm_2[bfactor.index(best_b)]:
                x_val = [-1* val for val in x_val]
                alt = True

            for x, b in zip(x_val, b_val):
                if x < 0.0:
                    alt_data.append(-1*x)
                    alt_b.append(b)
                else:
                    dom_data.append(x)
                    dom_b.append(b)

            dom_name = df.columns[4].split('_')[0]
            alt_name = df.columns[2].split('_')[0]


        number_line_dict["x_values"] = x_val
        number_line_dict["bfactor"] = b_val
        number_line_dict["name"] = name
        number_line_dict["y_values"] = [0] * len(number_line_dict["x_values"])
        df_number_line = pd.DataFrame.from_dict(number_line_dict)

        print(df_number_line.sort_values(by=['x_values'], ascending=True))

        df_number_line = df_number_line.sort_values(by=['bfactor'], ascending=True)
        im = axes[idx].scatter(df_number_line['x_values'], df_number_line['y_values'], c=df_number_line['bfactor'], vmin=0.5, vmax=0.9)

        axes[idx].set_xlim([-1.1, 1.1])

        # Build a rectangle in axes coords to align the text
        left, width = 0, 1
        bottom, height = 0, 1
        right = left + width
        top = bottom + height
        ax = plt.gca()
        p = plt.Rectangle((left, bottom), width, height, fill=False, alpha=0.0)
        p.set_transform(axes[idx].transAxes)
        p.set_clip_on(False)
        axes[idx].add_patch(p)

        if alt == False:
            axes[idx].text(left, 0.65 * (bottom + top), '{:}'.format(dom_name), verticalalignment='top', horizontalalignment='left', transform=axes[idx].transAxes, fontsize=12, rotation=0, fontweight='heavy')
            axes[idx].text(right, 0.65 * (bottom + top), '{:}'.format(alt_name), verticalalignment='top', horizontalalignment='right', transform=axes[idx].transAxes, fontsize=12, rotation=0, fontweight='heavy')
        else:
            axes[idx].text(left, 0.65 * (bottom + top), '{:}'.format(alt_name), verticalalignment='top', horizontalalignment='left', transform=axes[idx].transAxes, fontsize=12, rotation=0, fontweight='heavy')
            axes[idx].text(right, 0.65 * (bottom + top), '{:}'.format(dom_name), verticalalignment='top', horizontalalignment='right', transform=axes[idx].transAxes, fontsize=12, rotation=0, fontweight='heavy')


    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_ticks([0.5,0.7,0.9])
    cbar.set_ticklabels(["≤0.5","0.7","≥0.9"])

    # plt.show()
    plt.savefig("tm_all.png", dpi=300)
    plt.clf()
    plt.close()


    #create the heatmap
    dom_x = np.array(dom_data)
    dom_y = np.array(dom_b)

    print(len([i for i in range(len(dom_x)) if dom_x[i] >= 0.6 and dom_y[i] >=0.7])/float(len(dom_x)))
    
    alt_x = np.array(alt_data)
    alt_y = np.array(alt_b)

    print(len([i for i in range(len(alt_x)) if alt_x[i] >= 0.6 and alt_y[i] >=0.7])/float(len(alt_x)))
    
    heatmap_dom, x_edges_dom, y_edges_dom = np.histogram2d(dom_x, dom_y, bins = 100, range=[[0,1],[0,1]])
    heatmap_dom[heatmap_dom == 0.] = np.nan
    heatmap_alt, x_edges_alt, y_edges_alt = np.histogram2d(alt_x, alt_y, bins = 100, range=[[0,1],[0,1]])
    heatmap_alt[heatmap_alt == 0.] = np.nan
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))
    
    
    im1 = ax2.imshow(heatmap_dom.T, origin='lower', cmap='viridis', extent=[0,1,0,1], aspect='auto')
    ax1.set_title("Alternative Predictions")
    ax1.set_xlabel("TMscore")
    ax1.set_ylabel("plDDT")
    ax1.axvline(x=0.6, color='b', label='axvline - full height', c='black', ls=':', lw=1)
    ax1.axhline(y=0.7, color='b', label='axvline - full height', c='black', ls=':', lw=1)


    im2 = ax1.imshow(heatmap_alt.T, origin='lower', cmap='viridis', extent=[0,1,0,1], aspect='auto')
    ax2.set_title("Dominant Predictions")
    ax2.set_xlabel("TMscore")
    ax2.set_ylabel("plDDT")
    ax2.axvline(x=0.6, color='b', label='axvline - full height', c='black', ls=':', lw=1)
    ax2.axhline(y=0.7, color='b', label='axvline - full height', c='black', ls=':', lw=1)


    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im2, cax=cbar_ax, label='Count')

    plt.savefig("tm_plddt_hist_2d.png", dpi=600)
    plt.clf()
    plt.close()


def main():
    files = glob.glob("*.csv")
    files = sorted(files)

    for file in files:
        make_plot(file)
    combo_plot(files)

if __name__ == "__main__":
    main()
