<?php 
    require_once($path.'/include/libs/Filebase/src/Backup.php');
    require_once($path.'/include/libs/Filebase/src/Cache.php');
    require_once($path.'/include/libs/Filebase/src/Config.php');
    require_once($path.'/include/libs/Filebase/src/Database.php');
    require_once($path.'/include/libs/Filebase/src/Document.php');
    require_once($path.'/include/libs/Filebase/src/Filesystem.php');
    require_once($path.'/include/libs/Filebase/src/Predicate.php');
    require_once($path.'/include/libs/Filebase/src/QueryLogic.php');
    require_once($path.'/include/libs/Filebase/src/Query.php');
    require_once($path.'/include/libs/Filebase/src/SortLogic.php');
    require_once($path.'/include/libs/Filebase/src/Validate.php');
    require_once($path.'/include/libs/Filebase/src/Filesystem/FilesystemException.php');
    require_once($path.'/include/libs/Filebase/src/Filesystem/ReadingException.php');
    require_once($path.'/include/libs/Filebase/src/Filesystem/SavingException.php');
    require_once($path.'/include/libs/Filebase/src/Format/FormatException.php');
    require_once($path.'/include/libs/Filebase/src/Format/DecodingException.php');
    require_once($path.'/include/libs/Filebase/src/Format/EncodingException.php');
    require_once($path.'/include/libs/Filebase/src/Format/FormatInterface.php');
    require_once($path.'/include/libs/Filebase/src/Format/Json.php');
    require_once($path.'/include/libs/Filebase/src/Format/Yaml.php'); 

    // open the used databases
    $users = new \Filebase\Database([
        'dir' => $path.'/storage/users'
    ]);

    $categories = new \Filebase\Database([
        'dir' => $path.'/storage/categories'
    ]);

?>