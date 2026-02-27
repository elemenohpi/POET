"""Legacy settings module.

This file is no longer used. Configuration is now handled by
eletility.ConfigParser, which reads config.ini into a flat dict
passed through the application. Kept for reference only.
"""

# This module previously read from a `config` variable that was never
# defined at module level, making it broken on import. All config
# access now goes through the `config` dict created by cli.manage_input().
