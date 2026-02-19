import logging


def is_reduced(word):
    from schubmult import Permutation
    return Permutation.ref_product(*word).inv == len(word)

def little_bump(word, i, depth=0):
    from schubmult import Permutation
    indent = "  " * depth
    logging.debug(f"{indent}little_bump called with word={word}, i={i} (depth={depth})")
    
    if depth > 100:
        logging.error(f"Maximum recursion depth exceeded!")
        return None
    
    if i - 1 < 0 or i - 1 >= len(word):
        logging.info(f"{indent}❌ Index {i} out of bounds for word {word} (valid: 1-{len(word)}).")
        return None
    
    hat_word = word[:i - 1] + word[i:]
    logging.debug(f"{indent}hat_word = {hat_word}")
    
    if not is_reduced(hat_word):
        logging.info(f"{indent}❌ Hatted word {hat_word} is NOT reduced.")
        return None
    
    if word[i - 1] > 1:
        new_word = [*hat_word[:i - 1], word[i - 1] - 1, *hat_word[i - 1:]]
        logging.debug(f"{indent}Case 1: word[{i-1}]={word[i-1]} > 1, new_word = {new_word}")
    else:
        new_word = [*[a + 1 for a in hat_word[:i - 1]], word[i - 1], *[a + 1 for a in hat_word[i - 1:]]]
        logging.debug(f"{indent}Case 2: word[{i-1}]={word[i-1]} <= 1, new_word = {new_word}")
    
    if is_reduced(new_word):
        logging.debug(f"{indent}✓ new_word {new_word} is reduced, returning it")
        return new_word
    
    logging.debug(f"{indent}new_word {new_word} is NOT reduced, finding j...")
    
    # find j
    if is_reduced(new_word[:i]):
        logging.debug(f"{indent}new_word[:i] = {new_word[:i]} is reduced, searching forward")
        j = i + 1
        while j <= len(new_word) and is_reduced(new_word[:j]):
            logging.debug(f"{indent}  checking new_word[:{j}] = {new_word[:j]}")
            j += 1
        logging.info(f"{indent}Found j={j} (forward search), recursing with little_bump(new_word, {j + 1})")
        return little_bump(new_word, j + 1, depth + 1)
    
    logging.debug(f"{indent}new_word[:i] = {new_word[:i]} is NOT reduced, searching backward")
    j = i - 1
    while j >= 0 and is_reduced(new_word[j:]):
        logging.debug(f"{indent}  checking new_word[{j}:] = {new_word[j:]}")
        j -= 1
    logging.info(f"{indent}Found j={j} (backward search), recursing with little_bump(new_word, {j + 1})")
    return little_bump(new_word, j + 1, depth + 1)


if __name__ == "__main__":
    from schubmult import *
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(levelname)s: %(message)s'
    )
    while True:
        try:
            word = [int(x) for x in input("Enter a reduced word (space-separated): ").split()]
            print(f"Checking if word {word} is reduced...")
            if not is_reduced(word):
                print(f"WARNING: Input word {word} is NOT reduced!")
            else:
                print(f"Word {word} is reduced (length = {len(word)})")
            
            i = int(input(f"Enter an index between 1 and {len(word)}: "))
            print(f"\nCalling little_bump(word={word}, i={i})")
            print("="*60)
            
            new_word = little_bump(word, i)
            
            print("="*60)
            print("Result:")
            if new_word is None:
                print("❌ little_bump returned None - no little bump possible.")
            else:
                print(f"✓ New word: {new_word}")
                print(f"  Is reduced: {is_reduced(new_word)}")
            print()
        except (KeyboardInterrupt, EOFError):
            print("\nExiting...")
            break
        except Exception as e:
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()