#!/bin/sh

# Define configuration file names
CONFIG_ENV_FNAME="config_env.sh"
UNCONFIG_ENV_FNAME="unconfig_env.sh"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Source the paths.sh file
. "$SCRIPT_DIR/paths.sh"

# Generate config_env.sh
cat << EOF > $CONFIG_ENV_FNAME
#!/bin/sh

# Determine the project root directory based on the location of this script
LRGSG_ROOT="\$(cd "\$(dirname "\$0")" && pwd)"

# Export directories
export LRGSG_ROOT
EOF

for path in "${paths[@]}"; do
    echo "export ${path}" >> $CONFIG_ENV_FNAME
done

# Generate unconfig_env.sh
cat << EOF > $UNCONFIG_ENV_FNAME
#!/bin/sh

# Unset directories
EOF

for path in "${paths[@]}"; do
    var_name=$(echo "${path}" | cut -d '=' -f 1)
    echo "unset ${var_name}" >> $UNCONFIG_ENV_FNAME
done

echo "Scripts $CONFIG_ENV_FNAME and $UNCONFIG_ENV_FNAME have been generated."
