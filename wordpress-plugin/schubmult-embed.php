<?php
/**
 * Plugin Name: Schubmult Embed
 * Description: Provides a [schubmult] shortcode that embeds the schubmult web
 *              widget (hosted separately) via an iframe. Configure the URL
 *              under Settings → Schubmult Embed.
 * Version:     1.0.0
 * Author:      schubmult
 * License:     MIT
 */

if (!defined('ABSPATH')) { exit; }

const SCHUBMULT_OPTION = 'schubmult_embed_url';
const SCHUBMULT_DEFAULT_URL = 'https://example.pythonanywhere.com/embed';

/**
 * Shortcode: [schubmult height="640" width="100%"]
 *
 * Renders an <iframe> pointing at the configured embed URL.
 */
function schubmult_embed_shortcode($atts) {
    $atts = shortcode_atts(array(
        'height' => '640',
        'width'  => '100%',
        'flavor' => '',     // optional: py | double; appended as ?flavor=...
    ), $atts, 'schubmult');

    $url = trim((string) get_option(SCHUBMULT_OPTION, SCHUBMULT_DEFAULT_URL));
    if (empty($url)) {
        return '<em>Schubmult embed URL is not configured. See Settings → Schubmult Embed.</em>';
    }

    $src = esc_url($url);
    if (!empty($atts['flavor']) && in_array($atts['flavor'], array('py', 'double'), true)) {
        $sep = (strpos($src, '?') === false) ? '?' : '&';
        $src .= $sep . 'flavor=' . rawurlencode($atts['flavor']);
    }

    $height = preg_replace('/[^0-9]/', '', (string) $atts['height']);
    if ($height === '') { $height = '640'; }
    $width  = esc_attr($atts['width']);

    return sprintf(
        '<iframe src="%s" width="%s" height="%s" loading="lazy" '
        . 'style="border: 1px solid #ddd; border-radius: 6px;" '
        . 'sandbox="allow-scripts allow-same-origin allow-forms"></iframe>',
        $src, $width, $height
    );
}
add_shortcode('schubmult', 'schubmult_embed_shortcode');

/**
 * Settings page.
 */
function schubmult_embed_admin_menu() {
    add_options_page(
        'Schubmult Embed',
        'Schubmult Embed',
        'manage_options',
        'schubmult-embed',
        'schubmult_embed_settings_page'
    );
}
add_action('admin_menu', 'schubmult_embed_admin_menu');

function schubmult_embed_settings_init() {
    register_setting('schubmult_embed', SCHUBMULT_OPTION, array(
        'type'              => 'string',
        'sanitize_callback' => 'esc_url_raw',
        'default'           => SCHUBMULT_DEFAULT_URL,
    ));
}
add_action('admin_init', 'schubmult_embed_settings_init');

function schubmult_embed_settings_page() {
    if (!current_user_can('manage_options')) { return; }
    $url = (string) get_option(SCHUBMULT_OPTION, SCHUBMULT_DEFAULT_URL);
    ?>
    <div class="wrap">
        <h1>Schubmult Embed</h1>
        <p>Set the URL of your hosted schubmult Flask app's <code>/embed</code> endpoint.
           Then drop <code>[schubmult]</code> into any post or page.</p>
        <form action="options.php" method="post">
            <?php settings_fields('schubmult_embed'); ?>
            <table class="form-table" role="presentation">
                <tr>
                    <th scope="row"><label for="schubmult_embed_url">Embed URL</label></th>
                    <td>
                        <input id="schubmult_embed_url" name="<?php echo esc_attr(SCHUBMULT_OPTION); ?>"
                               type="url" value="<?php echo esc_attr($url); ?>"
                               class="regular-text" required>
                        <p class="description">e.g. <code>https://yourname.pythonanywhere.com/embed</code></p>
                    </td>
                </tr>
            </table>
            <?php submit_button(); ?>
        </form>
        <h2>Shortcode reference</h2>
        <pre>[schubmult]
[schubmult height="500"]
[schubmult flavor="double" height="700"]</pre>
        <p>The host serving the embed must include your WordPress site in its
           <code>SCHUBMULT_ALLOWED_ORIGINS</code> environment variable so the
           iframe is allowed.</p>
    </div>
    <?php
}
