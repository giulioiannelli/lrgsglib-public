#!/bin/sh

# Define configuration file names
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONFIG_ENV_FNAME="config_env.sh"
UNCONFIG_ENV_FNAME="unconfig_env.sh"
CONFIG_ENV_PATH="$SCRIPT_DIR/$CONFIG_ENV_FNAME"
UNCONFIG_ENV_PATH="$SCRIPT_DIR/$UNCONFIG_ENV_FNAME"
# Source the paths.sh file
. "$SCRIPT_DIR/paths.sh"

# Generate config_env.sh
cat << EOF > $CONFIG_ENV_PATH
#!/bin/sh

# Determine the project root directory based on the location of this script
LRGSG_ROOT="\$(cd "\$(dirname "\$0")" && pwd)"

# Export directories
export LRGSG_ROOT
EOF

for path in "${paths[@]}"; do
    echo "export ${path}" >> $CONFIG_ENV_PATH
done

# Generate unconfig_env.sh
cat << EOF > $UNCONFIG_ENV_PATH
#!/bin/sh

# Unset directories
EOF

for path in "${paths[@]}"; do
    var_name=$(echo "${path}" | cut -d '=' -f 1)
    echo "unset ${var_name}" >> $UNCONFIG_ENV_PATH
done

echo "Scripts $CONFIG_ENV_FNAME and $UNCONFIG_ENV_FNAME have been generated."
