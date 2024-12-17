import numpy as np
import matplotlib.pyplot as plt

def read_data(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
#            content = file.read()
#            tracker_only_phrase = "YES track ref and NO muon ref"
#            muon_only_phrase = "NO track ref and YES muon ref"
#            tracker_muon_phrase = "YES track ref and YES muon ref"
#            tracker_only_count = content.lower().count(tracker_only_phrase.lower())
#            muon_only_count = content.lower().count(muon_only_phrase.lower())
#            tracker_muon_count = content.lower().count(tracker_muon_phrase.lower())
#            print(f"Tracker only: {tracker_only_count}")
#            print(f"Muon only: {muon_only_count}")
#            print(f"Tracker + Muon: {tracker_muon_count}")
            
            lines = file.readlines()
            TrackMuonref = 0
            TrackMuonref_passcut = 0
            Trackref = 0
            Trackref_passcut = 0
            Muonref = 0
            Muonref_passcut = 0
            exit1 = 0
            exit2 = 0
            
            TrackMuon_track_pT_array = []
            TrackMuon_muon_pT_array = []
            Track_only_pT_array = []
            Muon_only_pT_array = []
            
            for line in lines:
                if line.startswith("Track + muon ref found"):
                    TrackMuonref += 1
                    TrackMuon_track_pT = line.split("Pt: ")[-2].strip(",  Muon ")
                    TrackMuon_muon_pT = line.split("Pt: ")[-1].strip("\n")
                    if float(TrackMuon_track_pT) > 45.0 and float(TrackMuon_muon_pT) > 5.0:
                        TrackMuonref_passcut += 1
                    #else:
                    #    print(f"Track pT: {TrackMuon_track_pT}, Muon pT: {TrackMuon_muon_pT}") 
                if line.startswith("Track ref only found"):
                    Trackref += 1
                    Track_only_pT = line.split("Pt: ")[-1].strip("\n")
                    Track_only_pT_array.append(Track_only_pT)
                    if float(Track_only_pT) > 45.0:
                        Trackref_passcut += 1
                    #else:
                    #    print(f"Track pT: {Track_only_pT}")
                if line.startswith("Muon ref only found"):
                    Muonref += 1
                    Muon_only_pT = line.split("Pt: ")[-1].strip("\n")
                    Muon_only_pT_array.append(Muon_only_pT)
                    if float(Muon_only_pT) > 5.0:
                        Muonref_passcut += 1
                    #else:
                    #    print(f"Muon pT: {Muon_only_pT}")
                if line.startswith("potential exit 1"):
                    exit1 += 1
                if line.startswith("potential exit 2"):
                    exit2 += 1
                total_exits = exit1 + exit2
                    
            print(f"{TrackMuonref_passcut}/{TrackMuonref} track + muon candidates pass cut")
            print(f"{Trackref_passcut}/{Trackref} track-only candidates pass cut")
            print(f"{Muonref_passcut}/{Muonref} muon-only candidates pass cut")
            print(f"Potential exit 1 for muon only references: {exit1}")
            print(f"Potential exit 2 for muon only references: {exit2}")
            print(f"Total exits (1 + 2) for muon only references: {total_exits}")
            
#            plt.hist(Track_only_pT_array, bins=50, alpha=0.7, color='blue', edgecolor='black')
#            plt.title('Track Only pT')
#            plt.xlabel('pT [GeV]')
#            plt.ylabel('')
#            #plt.grid(axis='y', alpha=0.75)
#            plt.show()
            
    except FileNotFoundError:
        print(f"The file '{file_path}' does not exist.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


#def plot_histogram(data):
#    plt.hist(data, 20, alpha=0.7, color='blue')
#    plt.title('Histogram')
#    plt.xlabel('muon pT')
#    #plt.ylabel('Frequency')
#    plt.grid(True)
#    plt.show()

def main():
    file_path = input("Enter the path to the .txt file: ")
 #   number_bins = input ("Enter the number of bins: ")
    data = read_data(file_path)
    
#    if data is not None:
#        plot_histogram(data, number_bins)

if __name__ == "__main__":
    main()
