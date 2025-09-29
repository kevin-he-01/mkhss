# Analyze passphrase files

import os
import math
import random
import tqdm
import pickle

from matplotlib import pyplot as plt

# Line Format:
# 12345<TAB>word

os.chdir(os.path.dirname(__file__))

def neighborhood_size(wordlist, word, threshold_hamming):
    # Calculate the neighborhood size of a wordlist file based on Hamming distance.
    """
    Calculate the neighborhood size of a wordlist file based on Hamming distance.
    Args:
        wordlist (list): A list of words from the wordlist file.
        word (str): The word to analyze.
        threshold_hamming (int): The Hamming distance threshold.
    Returns:
        int: The neighborhood size.
    """

    neighbor_count = 0
    for other_word in wordlist:
        assert len(other_word) == len(word), f'Word {other_word} has length {len(other_word)}, expected {len(word)}'
        distance = sum(c1 != c2 for c1, c2 in zip(word, other_word))
        # if word.startswith('recall'):
        #     print('E', word, other_word, distance)
        if distance <= threshold_hamming:
            neighbor_count += 1
    return neighbor_count

def max_neighborhood_arbitrary_word(wordlist: list[str], L: int, charset, threshold_hamming: int):
    # Calculate the maximum neighborhood size for an arbitrary word (not necessarily in the wordlist), where neighborhood size is defined as the number of words in the wordlist that are within a certain Hamming distance from the word.
    assert threshold_hamming == 2, 'This function is only implemented for threshold_hamming = 2'
    neighbor_counts: dict[str, int] = {}
    for word in tqdm.tqdm(wordlist):
        assert len(word) == L, f'Word {word} has length {len(word)}, expected {L}'
        neighbors: set[str] = set()
        for i in range(L): # Assumes threshold_hamming = 2
            for j in range(i + 1, L):
                new_word = list(word)
                for c_i in charset:
                    for c_j in charset:
                        # Create a new word by flipping the i-th and j-th characters
                        new_word[i] = c_i
                        new_word[j] = c_j
                        neighbors.add(''.join(new_word))
        for neighbor in neighbors:
            if neighbor not in neighbor_counts:
                neighbor_counts[neighbor] = 0
            neighbor_counts[neighbor] += 1
    # print(neighbor_counts)
    max_entry = max(neighbor_counts.items(), key=lambda x: x[1])
    print('Max entry: ', max_entry)
    return max_entry[1], neighbor_counts

def greedy_filter(wordlist: list[str], L: int, charset, threshold_hamming: int, nsize_max: int):
    assert threshold_hamming == 2, 'This function is only implemented for threshold_hamming = 2'
    try:
        with open('bad_neighbors.pkl', 'rb') as f:
            bad_neighbors = pickle.load(f)
    except FileNotFoundError:
        print('No bad neighbors file found, calculating bad neighbors...')
        bad_neighbors = []
        _, neighbor_counts = max_neighborhood_arbitrary_word(wordlist, L, charset, threshold_hamming)
        for neighbor, count in tqdm.tqdm(neighbor_counts.items()):
            if count > nsize_max:
                assert neighbor not in wordlist, f'Neighbor {neighbor} is in the wordlist, but has neighborhood size {count} > {nsize_max}'
                bad_neighbors.append((neighbor, count))
        print(f'Found {len(bad_neighbors)} bad neighbors with neighborhood size > {nsize_max}')
        # Pickle bad neighbors to a file
        with open('bad_neighbors.pkl', 'wb') as f:
            pickle.dump(bad_neighbors, f)
    print('Bad neighbors:')
    worst_case_num_words_removed = 0
    for neighbor, count in bad_neighbors:
        print(f'  {neighbor} ({count})')
        worst_case_num_words_removed += count - nsize_max
    print(f'Worst case number of words removed: {worst_case_num_words_removed} out of {len(wordlist)} words ({worst_case_num_words_removed / len(wordlist) * 100:.2f}%)')
    bad_neighbor_tracker: dict[str, int] = {neighbor: count for neighbor, count in bad_neighbors}
    print('Filtering wordlist...')
    removed_words = []
    # for word in tqdm.tqdm(wordlist):
    while bad_neighbor_tracker:
        # print('*', end='')
        offense_counts: list[tuple[str, int]] = []
        for word in wordlist:
            assert len(word) == L, f'Word {word} has length {len(word)}, expected {L}'
            if word in removed_words:
                # print(f'Skipping removed word {word}')
                continue
            offense_counts.append((word, neighborhood_size(list(bad_neighbor_tracker.keys()), word, threshold_hamming)))
        del word
        # print(offense_counts[:10])
        max_word, max_offense_count = max(offense_counts, key=lambda x: x[1])
        if max_offense_count != 0:
            print(f'Max offense count: {max_offense_count} for word {max_word}')
        removed_words.append(max_word)
        offense_count_check = 0
        to_be_removed = set()
        for bad_neighbor in bad_neighbor_tracker:
            assert len(bad_neighbor) == L, f'Bad neighbor {bad_neighbor} has length {len(bad_neighbor)}, expected {L}'
            assert len(max_word) == L, f'Word {max_word} has length {len(max_word)}, expected {L}'
            distance = sum(c1 != c2 for c1, c2 in zip(max_word, bad_neighbor))
            # print('A', max_word, bad_neighbor, distance)
            if distance <= threshold_hamming:
                offense_count_check += 1
                bad_neighbor_tracker[bad_neighbor] -= 1
                assert bad_neighbor_tracker[bad_neighbor] >= nsize_max
                if bad_neighbor_tracker[bad_neighbor] == nsize_max:
                    # del bad_neighbor_tracker[bad_neighbor]
                    to_be_removed.add(bad_neighbor)
        for bad_neighbor in to_be_removed:
            del bad_neighbor_tracker[bad_neighbor]
        assert offense_count_check == max_offense_count, f'Offense count check failed: {offense_count_check} != {max_offense_count}'

    print(f'Number of words removed: {len(removed_words)} out of {len(wordlist)} words ({len(removed_words) / len(wordlist) * 100:.2f}%)')
    filtered_wordlist = [word.strip('\0') for word in wordlist if word not in removed_words]
    print(f'Filtered wordlist size: {len(filtered_wordlist)}')
    return filtered_wordlist

def analyze_passphrase_file(file_path, neighborhood_threshold=2, plot: bool = False, analyze_neighborhood: bool = False):
    """
    Analyze a passphrase file to determine its length and the number of unique words.
    
    Args:
        file_path (str): The path to the passphrase file.
    
    Returns:
        tuple: A tuple containing the length of the passphrase and the number of unique words.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        words = list()
        for line in lines:
            if '\t' not in line:
                continue
            parts = line.strip().split('\t')
            # if len(parts) > 1:
            #     words.append(parts[1])
            words.append(parts[-1])
    
    print('File:', file_path)
    print('First and last 10 words:')
    print(words[:10])
    print(words[-10:])
    print('')

    max_length = 0
    min_length = 0
    total_length = 0
    charset = {'\0'}
    
    for word in words:
        assert '\0' not in word, 'Word contains null character'
        charset.update(set(word))
        length = len(word)
        total_length += length
        if length > max_length:
            max_length = length
        if min_length == 0 or length < min_length:
            min_length = length
    
    avg_length = total_length / len(words) if words else -1
    entropy = math.log2(len(words)) if words else 0
    
    print('Statistics:')
    print(f'  Total words: {len(words)}')
    print(f'  Max length: {max_length}')
    print(f'  Min length: {min_length}')
    print(f'  Avg length: {avg_length:.2f}')
    print(f'  Entropy: {entropy:.2f} bits')
    print(f'  Charset size (including NUL): {len(charset)}')
    print(f'  Number of bits to represent charset: {math.ceil(math.log2(len(charset)))} bits')
    print('')

    print('Charset: ', repr(''.join(sorted(charset))))

    # Pad with '\0' to maximum length
    words_null_padded = [word.ljust(max_length, '\0') for word in words]
    print('First and last 10 words (null-padded):')
    print(words_null_padded[:10])
    print(words_null_padded[-10:])
    print('')

    print(f'Neighborhood size analysis: threshold = {neighborhood_threshold}')
    for word in random.sample(words_null_padded, min(10, len(words_null_padded))):
    # for word in ['dase\x00\x00\x00\x00\x00']:
        # Calculate neighborhood size with a threshold of neighborhood_threshold
        threshold_hamming = neighborhood_threshold
        n_size = neighborhood_size(words_null_padded, word, threshold_hamming)
        print(f'  Neighborhood size for "{word}": {n_size}')
        remaining_entropy = math.log2(len(words_null_padded) / n_size) if n_size > 0 else -1
        # Attacker need ~2^remiaining_entropy guesses to find a neighboring word
        print(f'    Remaining entropy: {remaining_entropy:.2f} bits')
    print('')

    if analyze_neighborhood or plot:
        remaining_entropies = []
        nsize_dist = {}
        n_sizes = []
        for word in tqdm.tqdm(words_null_padded):
            threshold_hamming = neighborhood_threshold
            n_size = neighborhood_size(words_null_padded, word, threshold_hamming)
            n_sizes.append(n_size)
            if n_size not in nsize_dist:
                nsize_dist[n_size] = 0
            nsize_dist[n_size] += 1
            remaining_entropy = math.log2(len(words_null_padded) / n_size) if n_size > 0 else -1
            remaining_entropies.append(remaining_entropy)

        print('Neighborhood size distribution:')
        for n_size, count in sorted(nsize_dist.items()):
            print(f'  Neighborhood size {n_size}: {count} words ({count / len(words_null_padded) * 100:.2f}%)')
        print('')

        # Plot this distribution as a histogram
        if analyze_neighborhood:
            plt.hist(n_sizes, bins=len(nsize_dist), edgecolor='black')
            plt.title(f'Neighborhood Size Distribution: {file_path}, Threshold: {neighborhood_threshold}')
            plt.xlabel('Neighborhood Size')
            plt.ylabel('Frequency')
            plt.grid(axis='y', alpha=0.75)
            plt.show()

    # Plot a histogram of remaining entropy
    if plot:
        # remaining_entropies = []
        # for word in words_null_padded:
        #     threshold_hamming = neighborhood_threshold
        #     n_size = neighborhood_size(words_null_padded, word, threshold_hamming)
        #     remaining_entropy = math.log2(len(words_null_padded) / n_size) if n_size > 0 else -1
        #     remaining_entropies.append(remaining_entropy)
        # print(remaining_entropies[:10])
        # Print average remaining entropy, and standard deviation
        avg_remaining_entropy = sum(remaining_entropies) / len(remaining_entropies) if remaining_entropies else 0
        print(f'Average remaining entropy: {avg_remaining_entropy:.2f} bits')
        # Do percentiles of remaining entropies
        percentiles = list(range(0, 100, 5))  # 0th to 100th percentiles
        remaining_entropies_sorted = sorted(remaining_entropies)
        for p in percentiles:
            percentile_value = remaining_entropies_sorted[int(len(remaining_entropies) * p / 100)]
            print(f'{p}th percentile of remaining entropy: {percentile_value:.2f} bits')
        plt.hist(remaining_entropies, bins=30, edgecolor='black')
        plt.title(f'Filename: {file_path}, Threshold: {neighborhood_threshold}')
        plt.xlabel('Remaining Entropy (bits)')
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
        plt.show()
    
    return words_null_padded, charset

def filter_neighborhood_size(wordlist_null_padded, threshold_hamming, nsize_max):
    """
    Filter a wordlist to only include words that have a neighborhood size below a certain threshold.
    
    Args:
        wordlist (list): A list of words from the wordlist file.
        threshold_hamming (int): The Hamming distance threshold.
    
    Returns:
        list: A filtered list of words.
    """
    filtered_words = []
    for i, word in enumerate(wordlist_null_padded):
        n_size = neighborhood_size(wordlist_null_padded, word, threshold_hamming)
        if n_size <= nsize_max:
            filtered_words.append(word.strip('\0'))  # Remove null padding for filtered words
        if i % 50 == 0:
            print(f'Processed {i} words, filtered {len(filtered_words)} words with neighborhood size <= {nsize_max}, Ratio: {len(filtered_words) / (i + 1):.2f}', end='\r')
    return filtered_words

# wordlists = 'diceware.wordlist.asc  eff_large_wordlist.txt  eff_short_wordlist_1.txt  eff_short_wordlist_2_0.txt'.split()
# for wordlist in wordlists:
#     analyze_passphrase_file(wordlist, plot=False)
# analyze_passphrase_file('eff_short_wordlist_1.txt', neighborhood_threshold=2)
# analyze_passphrase_file('diceware.wordlist.asc', neighborhood_threshold=1)

# words, _ = analyze_passphrase_file('eff_large_wordlist.txt', neighborhood_threshold=2, plot=False)
# print(words[:10])
# filtered = filter_neighborhood_size(words, 2, 10)
# print(f'Filtered {len(filtered)} words from {len(words)} words with neighborhood size < 100. Ratio: {len(filtered) / len(words):.2f}')
# # Save the filtered words to a file
# with open('filtered_eff_large_wordlist.txt', 'w') as f:
#     for word in filtered:
#         f.write('\t' + word + '\n')

# wordlist, charset = analyze_passphrase_file('filtered_eff_large_wordlist.txt', neighborhood_threshold=2, plot=False)
# # Max entry:  ('remal\x00\x00\x00\x00', 15)
# # Neighbors: fetal legal petal rebel recall rehab relax remake remark remix remold rival roman royal rural
# # max_neighbor, _ = max_neighborhood_arbitrary_word(wordlist, len(wordlist[0]), charset, 2)
# # print('Max neighborhood size:', max_neighbor)
# final_wordlist = greedy_filter(wordlist, len(wordlist[0]), charset, 2, 10)
# # Save the final wordlist to a file
# with open('final_filtered_eff_large_wordlist.txt', 'w') as f:
#     for word in final_wordlist:
#         f.write('\t' + word + '\n')

# wordlist, charset = analyze_passphrase_file('eff_large_wordlist.txt', neighborhood_threshold=2, plot=False)
# # Max entry:  ('dase\x00\x00\x00\x00\x00', 59)
# # Neighbors: bash cage cake cane cape case cash dab dad dares darn dart dash data dawn dice dime dish disk dole dose dove doze dude duke dupe dusk dust easel fade fame game gas gave hash hate lake lash last late mace name nape race rage rake rare rash rise ruse sage sake same sash take task wake wasp wise
# # Remaining entropy: 7.04 bits
# # max_neighbor, _ = max_neighborhood_arbitrary_word(wordlist, len(wordlist[0]), charset, 2)

wordlist, charset = analyze_passphrase_file('final_filtered_eff_large_wordlist.txt', neighborhood_threshold=2, plot=False)
max_neighbor, _ = max_neighborhood_arbitrary_word(wordlist, len(wordlist[0]), charset, 2)
print('NUM_NEIGHBORS_MAX =', max_neighbor)

# analyze_passphrase_file('eff_large_wordlist.txt', neighborhood_threshold=2, analyze_neighborhood=True)
# analyze_passphrase_file('filtered_eff_large_wordlist.txt', neighborhood_threshold=2, analyze_neighborhood=True)
# analyze_passphrase_file('final_filtered_eff_large_wordlist.txt', neighborhood_threshold=2, analyze_neighborhood=True, plot=True)
