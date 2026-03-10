"""
RNA Tokenizer for Viral Sequences
=================================
Specialized tokenizer for viral RNA sequences with support for:
1. Nucleotide-level tokenization (A, U, G, C, N)
2. K-mer tokenization for capturing local sequence patterns
3. BPE tokenization for learned subword units
4. Special tokens for sequence boundaries and masking
"""

from typing import List, Dict, Optional, Tuple, Union
import re
from collections import Counter, defaultdict
import torch
from torch import Tensor


class RNATokenizer:
    """Tokenizer for RNA sequences with multiple encoding strategies."""
    
    def __init__(self, strategy: str = "nucleotide", k: int = 3, vocab_size: int = 4096):
        """
        Args:
            strategy: "nucleotide", "kmer", or "bpe"
            k: k-mer size for k-mer tokenization
            vocab_size: vocabulary size for BPE tokenization
        """
        self.strategy = strategy
        self.k = k
        self.vocab_size = vocab_size
        
        # Special tokens
        self.special_tokens = {
            '[PAD]': 0,
            '[UNK]': 1, 
            '[MASK]': 2,
            '[CLS]': 3,
            '[SEP]': 4,
        }
        
        self.token_to_id = self.special_tokens.copy()
        self.id_to_token = {v: k for k, v in self.special_tokens.items()}
        self.next_id = len(self.special_tokens)
        
        if strategy == "nucleotide":
            self._init_nucleotide_vocab()
        elif strategy == "kmer":
            self._init_kmer_vocab()
        # BPE vocab will be built from data
            
    def _init_nucleotide_vocab(self):
        """Initialize nucleotide-level vocabulary."""
        nucleotides = ['A', 'U', 'G', 'C', 'N']  # N for ambiguous nucleotides
        for nt in nucleotides:
            self.token_to_id[nt] = self.next_id
            self.id_to_token[self.next_id] = nt
            self.next_id += 1
            
    def _init_kmer_vocab(self):
        """Initialize k-mer vocabulary (will be populated during training)."""
        # K-mer vocab built dynamically from sequences
        pass
    
    def clean_sequence(self, sequence: str) -> str:
        """Clean and normalize RNA sequence."""
        # Convert to uppercase and replace T with U for RNA
        sequence = sequence.upper().replace('T', 'U')
        # Keep only valid nucleotides
        sequence = re.sub(r'[^AUGCN]', 'N', sequence)
        return sequence
    
    def tokenize_nucleotide(self, sequence: str) -> List[str]:
        """Tokenize at nucleotide level."""
        return list(self.clean_sequence(sequence))
    
    def tokenize_kmer(self, sequence: str) -> List[str]:
        """Tokenize using k-mers with sliding window."""
        sequence = self.clean_sequence(sequence)
        if len(sequence) < self.k:
            return [sequence]  # Return whole sequence if shorter than k
        
        kmers = []
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]
            kmers.append(kmer)
            
            # Add k-mer to vocabulary if not present
            if kmer not in self.token_to_id:
                self.token_to_id[kmer] = self.next_id
                self.id_to_token[self.next_id] = kmer
                self.next_id += 1
                
        return kmers
    
    def build_bpe_vocab(self, sequences: List[str], vocab_size: int = None):
        """Build BPE vocabulary from sequences."""
        if vocab_size is None:
            vocab_size = self.vocab_size
            
        # Start with character-level tokens
        char_freq = Counter()
        for seq in sequences:
            seq = self.clean_sequence(seq)
            char_freq.update(seq)
            
        # Add characters to vocabulary
        for char, _ in char_freq.most_common():
            if char not in self.token_to_id:
                self.token_to_id[char] = self.next_id
                self.id_to_token[self.next_id] = char
                self.next_id += 1
        
        # Build BPE merges
        word_freqs = defaultdict(int)
        for seq in sequences:
            seq = self.clean_sequence(seq)
            # Split into characters with end-of-word marker
            word = ' '.join(seq) + ' </w>'
            word_freqs[word] += 1
        
        # Iteratively merge most frequent pairs
        while len(self.token_to_id) < vocab_size:
            pairs = defaultdict(int)
            
            # Count all adjacent pairs
            for word, freq in word_freqs.items():
                symbols = word.split()
                for i in range(len(symbols) - 1):
                    pairs[(symbols[i], symbols[i + 1])] += freq
            
            if not pairs:
                break
                
            # Find most frequent pair
            best_pair = max(pairs, key=pairs.get)
            
            # Merge the pair
            new_token = ''.join(best_pair)
            if new_token not in self.token_to_id:
                self.token_to_id[new_token] = self.next_id
                self.id_to_token[self.next_id] = new_token
                self.next_id += 1
            
            # Update word frequencies with merged token
            new_word_freqs = {}
            for word, freq in word_freqs.items():
                new_word = word.replace(' '.join(best_pair), new_token)
                new_word_freqs[new_word] = freq
            word_freqs = new_word_freqs
    
    def tokenize(self, sequence: str) -> List[str]:
        """Tokenize sequence based on strategy."""
        if self.strategy == "nucleotide":
            return self.tokenize_nucleotide(sequence)
        elif self.strategy == "kmer":
            return self.tokenize_kmer(sequence)
        elif self.strategy == "bpe":
            return self.tokenize_bpe(sequence)
        else:
            raise ValueError(f"Unknown tokenization strategy: {self.strategy}")
    
    def tokenize_bpe(self, sequence: str) -> List[str]:
        """Tokenize using BPE (simplified version)."""
        sequence = self.clean_sequence(sequence)
        # For now, fall back to character-level if BPE not fully implemented
        tokens = []
        for char in sequence:
            if char in self.token_to_id:
                tokens.append(char)
            else:
                tokens.append('[UNK]')
        return tokens
    
    def encode(self, sequence: str, add_special_tokens: bool = True) -> List[int]:
        """Convert sequence to token IDs."""
        tokens = self.tokenize(sequence)
        
        if add_special_tokens:
            tokens = ['[CLS]'] + tokens + ['[SEP]']
        
        ids = []
        for token in tokens:
            if token in self.token_to_id:
                ids.append(self.token_to_id[token])
            else:
                ids.append(self.token_to_id['[UNK]'])
        
        return ids
    
    def decode(self, ids: List[int]) -> str:
        """Convert token IDs back to sequence."""
        tokens = []
        for id in ids:
            if id in self.id_to_token:
                token = self.id_to_token[id]
                if token not in self.special_tokens:
                    tokens.append(token)
        
        if self.strategy == "nucleotide":
            return ''.join(tokens)
        elif self.strategy == "kmer":
            # Reconstruct from k-mers (simplified - may have overlaps)
            if not tokens:
                return ""
            result = tokens[0]
            for kmer in tokens[1:]:
                # Add only the last character of each k-mer to avoid full repetition
                result += kmer[-1] if len(kmer) > 0 else ""
            return result
        else:
            return ''.join(tokens)
    
    def encode_batch(self, sequences: List[str], max_length: int = None, 
                    padding: bool = True) -> Tuple[Tensor, Tensor]:
        """Encode batch of sequences with padding."""
        encoded = [self.encode(seq) for seq in sequences]
        
        if max_length is None:
            max_length = max(len(seq) for seq in encoded)
        
        # Truncate if necessary
        encoded = [seq[:max_length] for seq in encoded]
        
        if padding:
            # Pad sequences
            padded = []
            attention_masks = []
            
            for seq in encoded:
                pad_length = max_length - len(seq)
                padded_seq = seq + [self.token_to_id['[PAD]']] * pad_length
                mask = [1] * len(seq) + [0] * pad_length
                
                padded.append(padded_seq)
                attention_masks.append(mask)
            
            return torch.tensor(padded), torch.tensor(attention_masks)
        else:
            return torch.tensor(encoded), None
    
    @property
    def vocab_size_actual(self) -> int:
        """Get actual vocabulary size."""
        return len(self.token_to_id)
    
    @property
    def mask_token_id(self) -> int:
        """Get mask token ID."""
        return self.token_to_id['[MASK]']
    
    @property
    def pad_token_id(self) -> int:
        """Get padding token ID."""
        return self.token_to_id['[PAD]']


def load_fasta_sequences(filepath: str, max_sequences: int = None) -> List[Tuple[str, str]]:
    """Load sequences from FASTA file."""
    sequences = []
    current_header = None
    current_sequence = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_sequence)))
                    if max_sequences and len(sequences) >= max_sequences:
                        break
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_sequence = []
            else:
                current_sequence.append(line)
        
        # Add last sequence
        if current_header is not None:
            sequences.append((current_header, ''.join(current_sequence)))
    
    return sequences


if __name__ == "__main__":
    # Test tokenizer
    tokenizer = RNATokenizer(strategy="nucleotide")
    
    test_seq = "AUGCUGAUCCGAUN"
    print(f"Original: {test_seq}")
    
    tokens = tokenizer.tokenize(test_seq)
    print(f"Tokens: {tokens}")
    
    ids = tokenizer.encode(test_seq)
    print(f"IDs: {ids}")
    
    decoded = tokenizer.decode(ids)
    print(f"Decoded: {decoded}")
    
    print(f"Vocab size: {tokenizer.vocab_size_actual}")
    print(f"Mask token ID: {tokenizer.mask_token_id}")


