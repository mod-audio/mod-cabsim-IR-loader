######################################
#
# cabsim-IR-loader
#
######################################

# where to find the source code - locally in this case
MOD_CABSIM_IR_LOADER_SITE_METHOD = local
MOD_CABSIM_IR_LOADER_SITE = $($(PKG)_PKGDIR)/

# even though this is a local build, we still need a version number
# bump this number if you need to force a rebuild
MOD_CABSIM_IR_LOADER_VERSION = 1

# dependencies (list of other buildroot packages, separated by space)
MOD_CABSIM_IR_LOADER_DEPENDENCIES = lv2

# LV2 bundles that this package generates (space separated list)
MOD_CABSIM_IR_LOADER_BUNDLES = cabsim-IR-loader.lv2

# call make with the current arguments and path. "$(@D)" is the build directory.
MOD_CABSIM_IR_LOADER_TARGET_MAKE = $(TARGET_MAKE_ENV) $(TARGET_CONFIGURE_OPTS) $(MAKE) -C $(@D)/source


# build command
define MOD_CABSIM_IR_LOADER_BUILD_CMDS
	$(MOD_CABSIM_IR_LOADER_TARGET_MAKE)
endef

# install command
define MOD_CABSIM_IR_LOADER_INSTALL_TARGET_CMDS
	$(MOD_CABSIM_IR_LOADER_TARGET_MAKE) install DESTDIR=$(TARGET_DIR)
endef


# import everything else from the buildroot generic package
$(eval $(generic-package))
