echo "Syncing remote folder"
rsync -avzhre ssh --include 'CHI_TECH/***' \
--include 'CHI_RESOURCES/***' \
--include 'Modules/***' \
--include 'CHI_TEST/***' \
--include 'CMakeLists.txt' \
--include 'configure.sh' \
--exclude '*' ./ janv4@$GPU_MACHINE://home/janv4/Desktop/ChiTech/synced --delete
