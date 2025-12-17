#!/usr/bin/env python3
"""
swif2_gen.py - Universal Generator Dispatcher

Routes to appropriate generator (RBHG, SPIZG, etc.) based on config file.

Usage:
    python swif2_gen.py <config_file>
    python swif2_gen.py RBHG.config
    python swif2_gen.py SPIZG.config
    
    python swif2_gen.py --generator RBHG --config myconfig.config  # Explicit

Detection logic:
    1. Check for GENERATOR_TYPE in config file
    2. Infer from config filename (RBHG.config → RBHG)
    3. Use --generator command line argument

Config file format:
    Each generator has its own config format, but should include:
        GENERATOR_TYPE    RBHG   # or SPIZG, etc.
    
    See GeneratorConfigExamples/ for templates.
"""

import sys
import os
import re
import argparse

# Add generators directory to path so we can import generator modules
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GENERATORS_DIR = os.path.join(SCRIPT_DIR, "generators")
sys.path.insert(0, GENERATORS_DIR)


def detect_generator_type(config_file, explicit_generator=None):
    """
    Detect which generator to use based on config file or explicit argument.
    
    Priority:
        1. Explicit --generator argument
        2. GENERATOR_TYPE in config file
        3. Infer from config filename (e.g., RBHG.config → RBHG)
    
    Args:
        config_file (str): Path to configuration file
        explicit_generator (str): Explicitly specified generator (from --generator)
    
    Returns:
        str: Generator type (e.g., 'RBHG', 'SPIZG')
    
    Raises:
        ValueError: If generator type cannot be determined
    """
    
    # Priority 1: Explicit argument
    if explicit_generator:
        print(f"Using explicitly specified generator: {explicit_generator}")
        return explicit_generator.upper()
    
    # Priority 2: Read GENERATOR_TYPE from config
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Strip inline comments
                if '#' in line:
                    line = line.split('#')[0].strip()
                
                # Look for GENERATOR_TYPE
                parts = line.split()
                if len(parts) >= 2 and parts[0].upper() == 'GENERATOR_TYPE':
                    generator = parts[1].upper()
                    print(f"Detected generator from config file: {generator}")
                    return generator
    
    # Priority 3: Infer from filename
    basename = os.path.basename(config_file)
    
    # Pattern: RBHG.config, RBHG_debug.config, SPIZG.config, etc.
    match = re.match(r'^([A-Z]+)', basename.upper())
    if match:
        generator = match.group(1)
        print(f"Inferred generator from filename '{basename}': {generator}")
        return generator
    
    # Could not determine
    raise ValueError(
        f"Could not determine generator type from '{config_file}'.\n"
        f"Options:\n"
        f"  1. Add 'GENERATOR_TYPE RBHG' (or SPIZG) to your config file\n"
        f"  2. Name your config like 'RBHG.config' or 'SPIZG.config'\n"
        f"  3. Use: python swif2_gen.py --generator RBHG {config_file}"
    )


def get_generator_module(generator_type, config_file_abs):
    """
    Import and return the appropriate generator module.
    
    Sets sys.argv before import so generators can read config file path during module initialization.
    
    Args:
        generator_type (str): Generator type (e.g., 'RBHG', 'SPIZG')
        config_file_abs (str): Absolute path to config file
    
    Returns:
        module: The generator's Python module
    
    Raises:
        ImportError: If generator module not found
    """
    generator_map = {
        'RBHG': 'RBHG.swif2_RBHG',
        'SPIZG': 'SPIZG.swif2_SPIZG',
        # Add future generators here:
        # 'GENX': 'GENX.swif2_GENX',
    }
    
    if generator_type not in generator_map:
        available = ', '.join(generator_map.keys())
        raise ImportError(
            f"Unknown generator type: {generator_type}\n"
            f"Available generators: {available}"
        )
    
    module_path = generator_map[generator_type]
    
    # Save original argv
    original_argv = sys.argv.copy()
    
    try:
        # Set sys.argv BEFORE importing so the generator's module-level code
        # can parse command line arguments correctly
        sys.argv = [f'swif2_{generator_type}.py', '--config', config_file_abs]
        
        # Import the generator module
        # e.g., "from RBHG import swif2_RBHG"
        parts = module_path.split('.')
        module = __import__(module_path, fromlist=[parts[-1]])
        print(f"Loaded generator module: {module_path}")
        return module
    except ImportError as e:
        raise ImportError(
            f"Could not import generator module '{module_path}'.\n"
            f"Make sure generators/{generator_type}/swif2_{generator_type}.py exists.\n"
            f"Error: {e}"
        )
    finally:
        # Restore original argv
        sys.argv = original_argv


def run_generator(generator_module, config_file_abs, args):
    """
    Execute the generator's main function.
    
    Args:
        generator_module: The imported generator module (already loaded with correct config)
        config_file_abs (str): Absolute path to config file (for reference)
        args: Parsed command-line arguments (for future use)
    
    Note:
        The generator has already been executed during import since most generators
        run their code at module level rather than having a main() function.
    """
    
    # Check if generator has a main() function (most don't, they run on import)
    if hasattr(generator_module, 'main'):
        print(f"\n{'='*70}")
        print(f"Running {generator_module.__name__}.main()")
        print(f"{'='*70}\n")
        generator_module.main()
    else:
        # Generator ran on import
        print(f"\n{'='*70}")
        print(f"Generator {generator_module.__name__} executed")
        print(f"{'='*70}\n")


def main():
    """Main entry point for swif2_gen.py"""
    
    parser = argparse.ArgumentParser(
        description='Universal generator dispatcher for SWIF2 workflow',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python swif2_gen.py RBHG.config
    python swif2_gen.py SPIZG.config
    python swif2_gen.py --generator RBHG myconfig.config
    python swif2_gen.py --help

Config file should contain GENERATOR_TYPE or follow naming convention.
See GeneratorConfigExamples/ for templates.
        """
    )
    
    parser.add_argument('config_file', 
                       help='Path to generator configuration file')
    
    parser.add_argument('--generator', '-g',
                       help='Explicitly specify generator type (RBHG, SPIZG, etc.)',
                       default=None)
    
    args = parser.parse_args()
    
    # Validate config file exists
    if not os.path.exists(args.config_file):
        print(f"ERROR: Config file not found: {args.config_file}")
        sys.exit(1)
    
    # Convert to absolute path once at the start
    config_file_abs = os.path.abspath(args.config_file)
    
    try:
        # Step 1: Detect generator type
        generator_type = detect_generator_type(config_file_abs, args.generator)
        
        # Step 2: Load generator module (sets sys.argv before import so config loads correctly)
        generator_module = get_generator_module(generator_type, config_file_abs)
        
        # Step 3: Run generator (if it has a main function; otherwise already ran on import)
        run_generator(generator_module, config_file_abs, args)
        
    except (ValueError, ImportError) as e:
        print(f"\nERROR: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\nUNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
